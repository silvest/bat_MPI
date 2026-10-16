// Proposal persistence: correlated Gaussian + fixed coordinate, checked analytically.
#include <BAT/BCModel.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCLog.h>
#include <TH1.h>
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

static void require(bool ok, const std::string& why) { if (!ok) throw std::runtime_error(why); }
class Model : public BCModel {
public:
    double offset = 0.;
    unsigned callbacks = 0, samples = 0;
    std::vector<double> sum;
    Model() : BCModel("proposal_test") {
        AddParameter("x",-10.,10.); AddParameter("y",-10.,10.); AddParameter("fixed",-1.,1.);
        GetParameter(0).SetPrior(std::make_shared<BCGaussianPrior>(0.3,2.));
        GetParameter(1).SetPrior(std::make_shared<BCGaussianPrior>(0.5,3.));
        GetParameter(2).SetPriorConstant(); GetParameter(2).Fix(0.4);
    }
    double LogLikelihood(const std::vector<double>& p) override {
        const double x=p[0]-1., y=p[1]+1., det=2.-0.65*0.65;
        return offset-0.5*(2.*x*x-1.3*x*y+y*y)/det;
    }
    void MCMCCurrentPointInterface(const std::vector<double>& p,int c,bool) override {
        require(p.size()==3 && c>=0 && unsigned(c)<GetNChains(),"bad callback"); ++callbacks;
    }
    void MCMCUserIterationInterface() override {
        if(sum.empty()) sum.assign(GetNChains()*2,0.);
        ++samples;
        for(unsigned c=0;c<GetNChains();++c) {
            const ChainState& s=fMCMCStates[c];
            require(s.parameters[2]==0.4,"fixed coordinate lost");
            require(std::abs(s.log_likelihood-LogLikelihood(s.parameters))<1e-10,"cached likelihood not re-evaluated");
            require(std::abs(s.log_prior-LogAPrioriProbability(s.parameters))<1e-10,"cached prior not re-evaluated");
            sum[c*2]+=s.parameters[0]; sum[c*2+1]+=s.parameters[1];
        }
    }
    std::vector<double> kernel() const {
        std::vector<double> out(1,GetProposalFunctionDof());
        for(unsigned c=0;c<GetNChains();++c) {
            out.insert(out.end(),fMCMCProposalFunctionScaleFactor[c].begin(),fMCMCProposalFunctionScaleFactor[c].end());
            if(GetProposeMultivariate()) for(unsigned i=0;i<2;++i) for(unsigned j=0;j<2;++j) {
                out.push_back(fMultivariateProposalFunctionCovariance[c][i][j]);
                out.push_back(fMultivariateProposalFunctionCholeskyDecomposition[c][i][j]);
            }
        }
        return out;
    }
    bool untouched() const { return fMCMCStates.empty(); }
    void checkCurrentDensities() {
        for (const auto& state : fMCMCStates) {
            require(std::abs(state.log_likelihood-LogLikelihood(state.parameters))<1e-10,"saved likelihood was reused");
            require(std::abs(state.log_prior-LogAPrioriProbability(state.parameters))<1e-10,"saved prior was reused");
        }
    }
    void checkMeans() const {
        // Analytic mean with likelihood precision divided by T; prior precision unchanged.
        const double det=2.-0.65*0.65;
        for(unsigned c=0;c<GetNChains();++c) {
            const double t=GetChainTemperature(c);
            const double a=2./det/t+0.25, b=-0.65/det/t, d=1./det/t+1./9.;
            const double hx=2.65/det/t+0.3/4., hy=-1.65/det/t+0.5/9.;
            const double ex=(d*hx-b*hy)/(a*d-b*b), ey=(a*hy-b*hx)/(a*d-b*b);
            require(std::abs(sum[c*2]/samples-ex)<0.14,"wrong x posterior mean");
            require(std::abs(sum[c*2+1]/samples-ey)<0.14,"wrong y posterior mean");
        }
    }
};
static void configure(Model& m,bool multivariate,bool pt) {
    m.SetNChains(2);
    if(pt) m.SetParallelTempering({1.,3.,9.},2,5);
    m.SetProposeMultivariate(multivariate);
    m.SetProposalFunctionDof(3.);
    m.SetNIterationsPreRunMin(2000); m.SetNIterationsPreRunMax(3000);
    m.SetNIterationsPreRunCheck(500); m.SetNIterationsPreRunFactorized(multivariate?1000:0);
    m.SetNIterationsRun(30000); m.SetInitialPositionScheme(BCEngineMCMC::kInitRandomPrior);
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv); int rank,size; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
    try {
        const std::string mode=argc>1?argv[1]:"multivariate";
        const bool multi=mode!="factorized", pt=mode=="pt";
        Model worker;
        BCLog::SetLogLevelScreen(BCLog::nothing);
        if(rank) {
            double recv[4],ll;
            for(;;) {
                MPI_Scatter(nullptr,4,MPI_DOUBLE,recv,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
                if(recv[0]==-1.) break;
                ll=recv[0]==1.?worker.LogEval(std::vector<double>(recv+1,recv+4)):-INFINITY;
                MPI_Gather(&ll,1,MPI_DOUBLE,nullptr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            }
        } else {
            const std::string file="proposal_"+mode+".txt";
            Model trained; configure(trained,multi,pt); trained.SetRandomSeed(76431);
            trained.SetProposalFunctionOutput(file);
            trained.MarginalizeAll(BCIntegrate::kMargMetropolis);
            const auto kernel=trained.kernel();
            // All compatibility failures must happen before chain initialization/collectives.
            for(unsigned bad=0;bad<7;++bad) {
                Model m; configure(m,multi,pt);
                if(bad==0) m.GetParameter(0).SetLimits(-9.,10.);
                if(bad==1) m.GetParameter(2).Fix(0.5);
                if(bad==2) m.SetNChains(trained.GetNChains()+1);
                if(bad==3) m.SetProposeMultivariate(!multi);
                if(bad==4) m.GetParameter(0).SetName("renamed");
                std::string input=file;
                if(bad==5) { input=file+".truncated"; std::ifstream src(file); std::string header; std::getline(src,header); std::ofstream dst(input); dst<<header<<'\n'; }
                if(bad==6) {
                    if(pt) m.SetParallelTempering({1.,4.,9.},2,5);
                    else input=file+".missing";
                }
                m.LoadProposalFunction(input,0);
                bool rejected=false;
                try { m.MetropolisPreRun(); } catch(const std::runtime_error&) { rejected=true; }
                require(rejected && m.untouched(),"incompatible proposal not rejected before initialization");
            }
            Model loaded; configure(loaded,multi,pt); loaded.SetRandomSeed(97531);
            loaded.LoadProposalFunction(file,500);
            // Save a second copy only after the discarded warm-up has succeeded.
            loaded.SetProposalFunctionOutput(file+".reused");
            loaded.MarginalizeAll(BCIntegrate::kMargMetropolis);
            require(loaded.GetProposalFunctionLoaded(),"proposal not loaded");
            require(loaded.kernel()==kernel,"proposal changed during reload or fixed warm-up/production");
            require(loaded.GetNIterationsConvergenceGlobal()==-1,"old convergence result inherited");
            require(loaded.callbacks==loaded.GetNChains()*30500*(multi?1:2),"adaptation ran or warm-up count incorrect");
            require(loaded.GetStatistics().n_samples==60000,"warm-up/hot samples entered posterior statistics");
            require(loaded.GetMarginalizedHistogram(0)->GetEntries()==60000,"warm-up entered histograms");
            loaded.checkMeans();
            // Zero warm-up is explicitly supported without adaptation.
            Model zero; configure(zero,multi,pt); zero.LoadProposalFunction(file,0); zero.SetNIterationsRun(2);
            zero.MarginalizeAll(BCIntegrate::kMargMetropolis);
            require(zero.callbacks==zero.GetNChains()*2*(multi?1:2),"zero warm-up did not skip adaptation");
            require(zero.kernel()==kernel,"proposal did not round-trip exactly");
            // A new target is always evaluated, even when structural metadata match.
            // Check serial initialization (before any MPI proposal), with zero warm-up.
            Model changed; configure(changed,multi,pt); changed.offset=17.;
            changed.GetParameter(0).SetPrior(std::make_shared<BCGaussianPrior>(0.8,2.));
            changed.LoadProposalFunction(file,0);
            changed.MetropolisPreRun();
            require(changed.GetProposalFunctionLoaded(),"new likelihood could not initialize");
            changed.checkCurrentDensities();
            std::vector<double> stop(size*4,0.); double recv[4];
            for(int p=1;p<size;++p) stop[p*4]=-1.;
            if(size>1) MPI_Scatter(stop.data(),4,MPI_DOUBLE,recv,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
            std::cout<<"PASS "<<mode<<": exact proposal round-trip, fixed parameters, rejected incompatible files, discarded warm-up, analytic posterior"<<std::endl;
        }
    } catch(const std::exception& e) { std::cerr<<"rank "<<rank<<": "<<e.what()<<std::endl; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Finalize();
}
