// Test warmup lifecycle through real MPI evaluations and cold ROOT output.
#include <BAT/BCModel.h>
#include <BAT/BCLog.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
static void check(bool value, const std::string& reason) { if (!value) throw std::runtime_error(reason); }
class Model : public BCModel {
public:
    unsigned settling_calls=0, production_calls=0;
    std::vector<double> frozen;
    Model() : BCModel("staged_test") {
        AddParameter("x",-5.,5.); AddParameter("y",-5.,5.); AddParameter("fixed",-1.,1.);
        SetPriorConstantAll(); GetParameter(2).Fix(.4);
        AddObservable("x_squared",0.,25.);
    }
    double LogLikelihood(const std::vector<double>& p) override {
        return -.5*(p[0]*p[0]+p[1]*p[1]-.6*p[0]*p[1]);
    }
    void CalculateObservables(const std::vector<double>& p) override { GetObservable(0)=p[0]*p[0]; }
    std::vector<double> kernel() const {
        std::vector<double> v{GetProposalFunctionDof(),double(fMultivariateCovarianceUpdates)};
        for(const auto& s:fMCMCProposalFunctionScaleFactor) v.insert(v.end(),s.begin(),s.end());
        for(const auto& matrix:fMultivariateProposalFunctionCovariance)
            for(int i=0;i<matrix.GetNrows();++i) for(int j=0;j<matrix.GetNcols();++j) v.push_back(matrix[i][j]);
        for(const auto& matrix:fMultivariateProposalFunctionCholeskyDecomposition)
            for(int i=0;i<matrix.GetNrows();++i) for(int j=0;j<matrix.GetNcols();++j) v.push_back(matrix[i][j]);
        return v;
    }
    void MCMCCurrentPointInterface(const std::vector<double>&, int, bool) override {
        if(GetPreRunStatus().stage=="settling") {
            check(UseMultivariateProposal()==GetProposeMultivariate(),"settling switched back to factorized proposals");
            if(frozen.empty()) frozen=kernel();
            check(frozen==kernel(),"proposal changed during settling"); ++settling_calls;
        } else if(fMCMCPhase==kMainRun) {
            if(frozen.empty()) frozen=kernel();
            check(frozen==kernel(),"proposal changed in production"); ++production_calls;
        }
    }
    bool untouched() const { return fMCMCStates.empty(); }
};
class RejectedModel : public Model {
public:
    double ProposalFunction(unsigned, unsigned) override { return 1000.; }
};
static void configure(Model& m,bool joint,bool pt,unsigned settling,bool smooth=false) {
    m.SetNChains(2); if(pt)m.SetParallelTempering({1.,2.,4.},2,2);
    m.SetProposeMultivariate(joint); m.SetProposalFunctionDof(0.);
    m.SetPreRunPolicy(BCEngineMCMC::kPreRunStaged,settling);
    if(smooth)m.SetScaleAdaptation(BCEngineMCMC::kScaleSmooth);
    m.SetNIterationsPreRunMax(101); m.SetNIterationsPreRunMin(1500);
    m.SetNIterationsPreRunCheck(20);m.SetNIterationsPreRunFactorized(joint?50:0);
    m.SetNIterationsRun(501);m.SetNLag(3);m.SetRandomSeed(37291);
    m.SetRValueParametersCriterion(.1); // Deliberately impossible; not an adaptation gate.
    m.SetInitialPositionScheme(BCEngineMCMC::kInitRandomPrior);
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    try {
        Model worker;BCLog::SetLogLevelScreen(BCLog::nothing);
        if(rank) {
            double recv[4],value;
            for(;;) {
                MPI_Scatter(nullptr,4,MPI_DOUBLE,recv,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
                if(recv[0]==-1.)break;
                value=recv[0]==1.?worker.LogEval(std::vector<double>(recv+1,recv+4)):-INFINITY;
                MPI_Gather(&value,1,MPI_DOUBLE,nullptr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            }
        } else {
            const std::string mode=argc>1?argv[1]:"factorized";
            const bool joint=mode!="factorized",pt=mode=="pt";
            const bool smooth=argc>2 && std::string(argv[2])=="smooth";
            const std::string prefix=std::string(smooth?"smooth_":"staged_")+mode+"_";
            Model trained;configure(trained,joint,pt,37,smooth);
            check(trained.GetNIterationsPreRun()==0,"prerun counted before start");
            if(pt)trained.SetPTDiagnosticsOutput(prefix,17,0);
            trained.SetPreRunStatusOutput(prefix+"warmup_status.csv");
            trained.SetProposalFunctionOutput(prefix+"proposal.txt");
            trained.WriteMarkovChain(prefix+"cold.root","RECREATE",true,false);
            trained.MarginalizeAll(BCIntegrate::kMargMetropolis);
            const auto status=trained.GetPreRunStatus();
            check(status.scale_adaptation==(smooth?BCEngineMCMC::kScaleSmooth:BCEngineMCMC::kScaleLegacy),"scale method not recorded");
            if(smooth)check(status.factorized_scale_updates>0 && (!joint || status.joint_scale_updates>0),"smooth learning-rate counters missing");
            check(status.adaptation_sweeps==101 && status.settling_sweeps==37,"warmup budgets not exact");
            check(trained.GetNIterationsPreRun()==138,"warmup count excluded settling");
            check(status.stage=="complete" && status.stop_reason=="adaptation_budget_exhausted","incorrect warmup status");
            check(status.adaptation_complete && status.settling_complete && !status.proposal_loaded,"missing warmup completion");
            check(!status.posterior_convergence_assessed && !status.legacy_rvalue_passed,"warmup asserted convergence");
            check(trained.GetNIterationsConvergenceGlobal()==-1,"staged warmup left a legacy convergence claim");
            const unsigned kernels=joint?1:2;
            check(trained.settling_calls==trained.GetNChains()*37*kernels,"settling count or proposal kernel incorrect");
            check(trained.production_calls==trained.GetNChains()*501*kernels,"production has warmup callbacks or missing draws");
            check(trained.GetStatistics().n_samples==2*167,"warmup/lag/hot samples contaminated production statistics");
            check(trained.GetMarginalizedHistogram(0)->GetEntries()==2*167,"warmup contaminated histogram");
            trained.CloseOutputFile();
            TFile file((prefix+"cold.root").c_str(),"READ");
            auto tree=dynamic_cast<TTree*>(file.Get("staged_test_mcmc"));
            check(tree && tree->GetEntries()==2*167,"warmup/hot samples saved to cold ROOT tree");
            int phase;unsigned chain;tree->SetBranchAddress("Phase",&phase);tree->SetBranchAddress("Chain",&chain);
            for(Long64_t i=0;i<tree->GetEntries();++i) {tree->GetEntry(i);check(phase==1 && chain<2,"non-production or hot ROOT row");}
            file.Close();
            if(pt) {
                check(trained.GetSwapAttempts()[0]==250 && trained.GetSwapAttempts()[1]==250,"warmup swaps leaked into production");
                std::ifstream input(prefix+"pt_checkpoints.csv");std::string s((std::istreambuf_iterator<char>(input)),{});
                check(s.find("prerun,")!=std::string::npos && s.find("warmup,")!=std::string::npos && s.find("production,")!=std::string::npos,"staged diagnostics lost a phase");
            }
            Model loaded;configure(loaded,joint,pt,37);
            loaded.LoadProposalFunction(prefix+"proposal.txt",7);
            loaded.MarginalizeAll(BCIntegrate::kMargMetropolis);
            check(loaded.GetPreRunStatus().proposal_loaded && loaded.GetPreRunStatus().adaptation_sweeps==0 &&
                  loaded.GetPreRunStatus().settling_sweeps==7,"loaded proposal repeated staged adaptation/settling");
            check(loaded.GetPreRunStatus().stop_reason=="loaded_proposal","loaded proposal status incorrect");
            check(loaded.settling_calls==loaded.GetNChains()*7*kernels,"loaded proposal warmup count incorrect");
            check(loaded.kernel()==trained.kernel(),"saved proposal changed on load");
            Model zero;configure(zero,joint,pt,0,smooth);zero.MetropolisPreRun();
            check(zero.GetPreRunStatus().settling_complete && zero.settling_calls==0 && zero.GetNIterationsPreRun()==101,"zero settling ignored");
            // Configuration errors must fail before MPI collectives or initialization.
            for(unsigned problem=0;problem<4;++problem) {
                Model bad;configure(bad,true,pt,0);
                if(problem==0)bad.SetNIterationsPreRunMax(0);
                if(problem==1)bad.SetNIterationsPreRunCheck(0);
                if(problem==2)bad.SetNIterationsPreRunFactorized(101);
                if(problem==3)bad.SetPreRunPolicy(BCEngineMCMC::kPreRunStaged,std::numeric_limits<int>::max());
                bool rejected=false;try{bad.MetropolisPreRun();}catch(const std::invalid_argument&){rejected=true;}
                check(rejected && bad.untouched(),"invalid staged configuration reached initialization");
            }
            RejectedModel rejected;configure(rejected,false,pt,13);
            rejected.MetropolisPreRun();
            check(rejected.GetPreRunStatus().adaptation_complete && rejected.GetPreRunStatus().settling_complete,
                  "stuck proposals made staged warmup endless");
            check(rejected.GetPreRunStatus().adaptation_zero_acceptance==rejected.GetNChains()*2 &&
                  rejected.GetPreRunStatus().settling_zero_acceptance==rejected.GetNChains()*2,
                  "stuck proposal problems were hidden");
            check(!rejected.GetPreRunStatus().posterior_convergence_assessed,"stuck proposal reported convergence");
            bool lateRejected=false;
            try {trained.SetPreRunPolicy(BCEngineMCMC::kPreRunStaged,0);} catch(const std::logic_error&){lateRejected=true;}
            check(lateRejected,"warmup policy changed on an initialized engine");
            Model legacy;check(legacy.GetPreRunPolicy()==BCEngineMCMC::kPreRunLegacy,"legacy default changed");
            check(legacy.GetScaleAdaptation()==BCEngineMCMC::kScaleLegacy,"scale-adaptation default changed");
            legacy.SetNChains(2);legacy.SetMinimumEfficiency(0.);legacy.SetMaximumEfficiency(1.);
            legacy.SetRValueParametersCriterion(100.);legacy.SetNIterationsPreRunMin(20);
            legacy.SetNIterationsPreRunMax(1000);legacy.SetNIterationsPreRunCheck(10);
            legacy.SetProposeMultivariate(false);legacy.SetRandomSeed(4178);legacy.MetropolisPreRun();
            check(legacy.GetPreRunStatus().adaptation_sweeps<1000 && legacy.GetPreRunStatus().stop_reason=="legacy_criteria_met",
                  "legacy early stopping changed");
            Model copy(trained);check(copy.GetPreRunPolicy()==BCEngineMCMC::kPreRunStaged && copy.GetNIterationsSettling()==37,"copy lost staged settings");
            std::vector<double> stop(size*4,0.);double recv[4];for(int i=1;i<size;++i)stop[i*4]=-1.;
            if(size>1)MPI_Scatter(stop.data(),4,MPI_DOUBLE,recv,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
            std::cout<<"PASS "<<(smooth?"smooth ":"")<<mode<<": exact budgets, frozen kernels, discarded settling, cached proposals, legacy stopping, cold-only production\n";
        }
    } catch(const std::exception& e) { std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Finalize();
}
