// End-to-end MPI tests against numerical quadrature, including nonconstant priors.
#include <BAT/BCModel.h>
#include <BAT/BCDisjointFlatPrior.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCLog.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TRandom3.h>
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

static void check(bool ok, const std::string& message)
{ if (!ok) throw std::runtime_error(message); }

using Intervals = std::vector<std::pair<double, double>>;
static const Intervals support = {{1., 3.}, {-4., -3.}};

static bool inSupport(double x)
{ return (-4. <= x && x <= -3.) || (1. <= x && x <= 3.); }

static void checkDisjointPrior()
{
    BCDisjointFlatPrior prior(support);
    for (double x : {-4., -3.5, -3., 1., 2., 3.})
        check(std::abs(prior.GetPrior(x) - 1./3.) < 1.e-14, "wrong uniform density on support");
    for (double x : {-5., 0., 4.})
        check(prior.GetLogPrior(x) == -INFINITY, "finite prior outside support");
    check(std::abs(prior.GetIntegral() - 1.) < 1.e-14, "flat union prior is not normalized");
    check(std::abs(prior.GetIntegral(-3.5, 2.5) - 2./3.) < 1.e-14, "wrong clipped integral");
    check(prior.GetIntegral(-2., 0.) == 0., "nonzero mass in gap");
    check(std::abs(prior.GetMean() - 1./6.) < 1.e-14, "wrong length-weighted mean");

    const double width = std::ldexp(1., -20);
    BCDisjointFlatPrior narrow({{.25, .25 + width}, {3.75, 3.75 + 2 * width}});
    const double edges[] = {0., 1., 2., 4.};
    TH1D bins("disjoint_bins", "", 3, edges);
    BCH1D plot = narrow.GetBCH1D(&bins);
    TH1* h = plot.GetHistogram();
    check(h && std::abs(h->Integral(1, 1, "width") - 1./3.) < 1.e-14
            && h->GetBinContent(2) == 0. && std::abs(h->Integral(3, 3, "width") - 2./3.) < 1.e-14,
          "prior histogram lost mass from narrow intervals");

    TRandom3 random(13421);
    unsigned left = 0;
    for (unsigned i = 0; i < 10000; ++i) {
        const double x = prior.GetRandomValue(-3.5, 2.5, &random);
        check(inSupport(x) && -3.5 <= x && x <= 2.5, "invalid draw from clipped support");
        if (x < 0.) ++left;
    }
    check(std::abs(left / 10000. - 0.25) < 0.025, "wrong sampling weights on clipped intervals");

    for (const auto& invalid : std::vector<Intervals>{
             {}, {{0., 0.}, {1., 2.}}, {{2., 1.}, {3., 4.}},
             {{0., 2.}, {1., 3.}}, {{0., 1.}, {1., 3.}}, {{0., 1.}, {2., INFINITY}}}) {
        bool rejected = false;
        try { BCDisjointFlatPrior bad(invalid); }
        catch (const std::invalid_argument&) { rejected = true; }
        check(rejected, "invalid interval union accepted");
    }
}

class Mixture : public BCModel {
public:
    std::vector<double> sum, positive;
    unsigned samples = 0, callbacks = 0, main_start_callbacks = 0;
    unsigned gapProposals = 0;
    unsigned long long evaluations = 0;
    bool bimodal, disjoint;
    Mixture(bool multimodal = true, bool flatUnion = false)
        : BCModel("pt_test"), bimodal(multimodal), disjoint(flatUnion) {
        AddParameter("x", -10., 10.);
        if (disjoint)
            GetParameter(0).SetPrior(std::make_shared<BCDisjointFlatPrior>(support));
        else
            GetParameter(0).SetPrior(std::make_shared<BCGaussianPrior>(0.7, 3.));
    }
    double LogLikelihood(const std::vector<double>& p) override {
        if (disjoint) {
            check(inSupport(p[0]), "evaluated likelihood outside disjoint prior support");
            return 0.;
        }
        if (!bimodal) return -0.5 * std::pow((p[0]-1.)/0.8, 2);
        double a = std::log(0.3 / 0.45) - 0.5 * std::pow((p[0]+3.)/0.45, 2);
        double b = std::log(0.7 / 0.8) - 0.5 * std::pow((p[0]-3.)/0.8, 2);
        double m = std::max(a,b);
        return m + std::log(std::exp(a-m) + std::exp(b-m));
    }
    double LogEval(const std::vector<double>& p) override {
        ++evaluations;
        return BCModel::LogEval(p);
    }
    std::vector<PTDiagnosticPartition> DefinePTDiagnosticPartitions() const override {
        return {{"sign", {"negative", "positive"}}};
    }
    void ClassifyPTDiagnosticRegions(const std::vector<double>& p, std::vector<unsigned>& labels) const override {
        labels.assign(1, p.at(0)<0. ? 0 : 1);
    }
    void MCMCUserIterationInterface() override {
        if (sum.empty()) { sum.assign(GetNChains(),0.); positive.assign(GetNChains(),0.); }
        if (!samples) main_start_callbacks = callbacks - GetNChains();
        ++samples;
        for (unsigned c=0; c<GetNChains(); ++c) {
            const ChainState& s = fMCMCStates[c];
            if (disjoint) check(inSupport(s.parameters[0]), "stored sample outside support");
            sum[c] += s.parameters[0];
            positive[c] += s.parameters[0]>0.;
            check(std::abs(s.log_probability - s.log_prior - s.log_likelihood)<1e-10,
                  "cached probability fields are inconsistent after exchange");
        }
    }
    void MCMCCurrentPointInterface(const std::vector<double>& p, int c, bool accepted) override {
        check(c>=0 && unsigned(c)<GetNChains() && p.size()==1, "proposal callback out of bounds");
        if (disjoint && -3. < p[0] && p[0] < 1.) {
            ++gapProposals;
            check(!accepted, "accepted proposal in the gap");
        }
        ++callbacks;
    }
    std::pair<double,double> reference(double temperature) {
        if (disjoint) return {1./6., 2./3.}; // Flat likelihood, one unit of width left and two right.
        double norm=0., mean=0., pos=0.;
        const double dx=0.0005;
        for (double x=-10.+dx/2.; x<10.; x+=dx) {
            const std::vector<double> p(1,x);
            double w=std::exp(LogLikelihood(p)/temperature + LogAPrioriProbability(p));
            norm+=w; mean+=w*x; if(x>0.) pos+=w;
        }
        return std::make_pair(mean/norm,pos/norm);
    }
};

int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);
    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
    try {
        const std::string scenario=argc>1?argv[1]:"serial";
        const bool disabled=scenario=="disabled", two=scenario=="two_temperatures";
        const bool disjoint = scenario == "disjoint";
        Mixture m(!disabled, disjoint);
        BCLog::SetLogLevelScreen(BCLog::nothing);
        if(rank) {
            double recv[2], value;
            while(true) {
                MPI_Scatter(nullptr,2,MPI_DOUBLE,recv,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
                if(recv[0]==-1.) break;
                value=recv[0]==1.?m.LogEval(std::vector<double>(1,recv[1])):-INFINITY;
                MPI_Gather(&value,1,MPI_DOUBLE,nullptr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            }
        } else {
            if (disjoint) checkDisjointPrior();
            // Configuration validation must fail without starting any collectives.
            for(const auto& invalid : std::vector<std::vector<double>>{{},{2.,4.},{1.,1.},{1.,INFINITY},{1.,0.5}}) {
                bool failed=false;
                try { m.SetParallelTempering(invalid,2,5); } catch(const std::invalid_argument&) {failed=true;}
                check(failed,"invalid ladder accepted");
            }
            m.SetNChains(2);
            if(!disabled) m.SetParallelTempering(two?std::vector<double>{1.,9.}:std::vector<double>{1.,3.,9.,27.},2,5);
            if(scenario=="partial") m.SetPTDiagnosticsOutput("partial_stream_",1000,100);
            m.SetRandomSeed(71357);
            m.SetProposeMultivariate(scenario!="serial" && !disabled);
            m.SetProposalFunctionDof(0.);
            m.SetNIterationsPreRunMin(4000); m.SetNIterationsPreRunMax(5000);
            m.SetNIterationsPreRunCheck(500);
            if (scenario=="factorized_warmup") m.SetNIterationsPreRunFactorized(1000);
            m.SetNIterationsRun(40000);
            m.SetInitialPositionScheme(BCEngineMCMC::kInitRandomPrior);
            // Without a jump, all replicas would remain in the left interval.
            if (disjoint)
                m.SetInitialPositions(std::vector<std::vector<double>>(m.GetNChains(), {-3.5}));
            const std::string filename="pt_"+scenario+".root";
            m.WriteMarkovChain(filename,"RECREATE",true,false);
            m.MarginalizeAll(BCIntegrate::kMargMetropolis);
            check(m.GetStatistics().n_samples==80000,"hot samples contaminated posterior statistics");
            check(m.GetMarginalizedHistogram(0)->GetEntries()==80000,"hot samples contaminated histogram");
            if (disjoint) {
                check(m.gapProposals > 0, "no proposals tested the excluded gap");
                for (unsigned c = 0; c < m.GetNColdChains(); ++c)
                    check(std::abs(m.positive[c] / m.samples - 2./3.) < .05,
                          "cold chain failed to explore both intervals in proportion to their widths");
            }
            // Check every temperature, not only the posterior: this detects accidental prior tempering.
            for(unsigned t=0;t<(disabled?1:m.GetTemperatures().size());++t) {
                const double temp=disabled?1.:m.GetTemperatures()[t];
                const auto expected=m.reference(temp);
                const unsigned n=m.GetNColdChains(); double mean=0., pos=0.;
                for(unsigned l=0;l<n;++l) { mean+=m.sum[t*n+l]; pos+=m.positive[t*n+l]; }
                mean/=n*m.samples; pos/=n*m.samples;
                std::cout<<scenario<<" T="<<temp<<" mean="<<mean<<" expected="<<expected.first
                         <<" mass(x>0)="<<pos<<" expected="<<expected.second<<std::endl;
                check(std::abs(mean-expected.first)<0.20,"incorrect tempered mean");
                check(std::abs(pos-expected.second)<0.04,"incorrect mode weights");
            }
            for(unsigned t=0;t<m.GetSwapAttempts().size();++t) {
                const unsigned expected=(two?8000:4000)*2;
                check(m.GetSwapAttempts()[t]==expected,"incorrect exchange schedule or warm-up counters");
                check(m.GetSwapAccepts()[t]>0,"no accepted swaps");
            }
            check(m.callbacks-m.main_start_callbacks==m.GetNChains()*40000,"dropped or duplicate proposal callbacks");
            m.CloseOutputFile();
            TFile f(filename.c_str(),"READ");
            auto tree=dynamic_cast<TTree*>(f.Get("pt_test_mcmc"));
            check(tree && tree->GetEntries()==80000,"wrong number of saved cold samples");
            unsigned chain, replica; double x, ll, lp, logp;
            tree->SetBranchAddress("Chain",&chain); tree->SetBranchAddress("Replica",&replica); tree->SetBranchAddress("x",&x);
            tree->SetBranchAddress("LogLikelihood",&ll); tree->SetBranchAddress("LogPrior",&lp);
            tree->SetBranchAddress("LogProbability",&logp);
            for(Long64_t row=0;row<tree->GetEntries();++row) {
                tree->GetEntry(row);
                check(chain<2,"hot chain in output tree");
                check(replica<m.GetNChains() && replica%2==chain,"saved replica ID inconsistent with independent ladder");
                check(std::abs(logp-ll-lp)<1e-10,"wrong saved probability");
                check(std::abs(ll-m.LogLikelihood(std::vector<double>(1,x)))<1e-10,"likelihood cache mismatches saved state");
            }
            auto parameters=dynamic_cast<TTree*>(f.Get("pt_test_parameters"));
            unsigned saved_count=0;
            parameters->SetBranchAddress("nchains",&saved_count); parameters->GetEntry(0);
            check(saved_count==2,"parameter tree must describe cold chains only");
            f.Close();
            // Cold output can be reloaded as an ordinary posterior sample.
            Mixture loaded(!disabled, disjoint); loaded.LoadMCMC(filename,"","",false); loaded.Remarginalize(false);
            check(loaded.GetNChains()==2 && loaded.GetStatistics().n_samples==80000,"cold-chain reload failed");
            std::vector<double> stop(size*2,0.); double dummy[2];
            for(int p=1;p<size;++p) stop[2*p]=-1.;
            if(size>1) MPI_Scatter(stop.data(),2,MPI_DOUBLE,dummy,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
            std::cout<<"PASS "<<scenario<<std::endl;
        }
        std::vector<unsigned long long> evaluations(size);
        MPI_Gather(&m.evaluations, 1, MPI_UNSIGNED_LONG_LONG, evaluations.data(), 1,
                   MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        if (!rank) {
            for (unsigned r=0;r<evaluations.size();++r)
                check(r<std::min<unsigned>(size,m.GetNChains()) ? evaluations[r]>0 : evaluations[r]==0,
                      "likelihood dispatch failed to use every available rank up to the replica count");
            std::cout << "All available replica-evaluation ranks were exercised (including rank 0)\n";
        }
    } catch(const std::exception& e) {
        std::cerr<<"rank "<<rank<<": "<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    MPI_Finalize();
}
