// Analytic target, MPI lifecycle, bounded schedule updates and cache restoration.
#include <BAT/BCModel.h>
#include <BAT/BCLog.h>
#include "../src/PTTemperatureSchedule.h"
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void check(bool ok, const std::string& why) { if (!ok) throw std::runtime_error(why); }
class Model : public BCModel {
public:
    std::vector<double> frozen, production_temperatures;
    unsigned long long draws = 0;
    double sum = 0., sum2 = 0.;
    Model() : BCModel("temperature_test") {
        AddParameter("x", -20., 20.); AddParameter("y", -20., 20.);
        SetPriorGauss("x", 0., 2.); SetPriorGauss("y", 0., 2.);
    }
    double LogLikelihood(const std::vector<double>& p) override {
        return -2.*((p[0]-1.)*(p[0]-1.)+p[1]*p[1]);
    }
    std::vector<PTDiagnosticPartition> DefinePTDiagnosticPartitions() const override {
        return {PTDiagnosticPartition("sign", {"negative", "positive"})};
    }
    void ClassifyPTDiagnosticRegions(const std::vector<double>& p, std::vector<unsigned>& labels) const override {
        labels.assign(1,p[0]>=0. ? 1 : 0);
    }
    std::vector<double> kernel() const {
        std::vector<double> v{double(fMultivariateCovarianceUpdates)};
        for (const auto& scales : fMCMCProposalFunctionScaleFactor) v.insert(v.end(),scales.begin(),scales.end());
        for (const auto& matrix : fMultivariateProposalFunctionCholeskyDecomposition)
            for (int i=0;i<matrix.GetNrows();++i) for (int j=0;j<matrix.GetNcols();++j) v.push_back(matrix[i][j]);
        return v;
    }
    void MCMCCurrentPointInterface(const std::vector<double>&, int chain, bool) override {
        if (GetPreRunStatus().stage=="settling" || fMCMCPhase==kMainRun) {
            if (frozen.empty()) { frozen=kernel(); production_temperatures=GetTemperatures(); }
            check(frozen==kernel(), "proposals changed after adaptation");
            check(GetTemperatures()==production_temperatures, "temperatures changed after adaptation");
        }
        if (fMCMCPhase==kMainRun && chain<int(GetNColdChains())) {
            // The callback argument is the candidate, including rejected moves.
            const double x=fMCMCStates[chain].parameters[0];
            ++draws; sum+=x; sum2+=x*x;
        }
    }
    bool untouched() const { return fMCMCStates.empty(); }
};
static void configure(Model& model, bool joint) {
    model.SetParallelTempering({1.,1.01,2.,16.},2,2);
    model.SetPreRunPolicy(BCEngineMCMC::kPreRunStaged,31);
    model.SetScaleAdaptation(BCEngineMCMC::kScaleSmooth);
    model.SetPTTemperatureAdaptation(true,50,3,.5);
    model.SetProposeMultivariate(joint); model.SetNIterationsPreRunFactorized(80);
    model.SetNIterationsPreRunCheck(20); model.SetNIterationsPreRunMax(600);
    model.SetNIterationsRun(10001); model.SetNLag(2); model.SetRandomSeed(8328);
    model.SetProposalFunctionDof(0.); model.SetInitialPositionScheme(BCEngineMCMC::kInitRandomPrior);
}
static void schedule_tests() {
    const std::vector<double> t{1.,2.,4.};
    auto result=BATPT::UpdateSchedule(t,{.8,.2},.5);
    check(std::abs(result.temperatures[1]-1./.59375)<1.e-13,"barrier inversion/orientation incorrect");
    for (double r : {0.,.1,1.}) {
        result=BATPT::UpdateSchedule(t,{r,r},1.);
        check(result.temperatures==t,"uniform barrier changed the grid");
    }
    auto current=std::vector<double>{1.,1.00001,2.,3.,1000.};
    for (unsigned k=0;k<50;++k) {
        result=BATPT::UpdateSchedule(current,k%2 ? std::vector<double>{0.,0.,1.,0.} : std::vector<double>{1.,0.,0.,0.},1.);
        check(result.temperatures.front()==1. && result.temperatures.back()==1000.,"endpoint changed");
        for (unsigned i=1;i<current.size();++i) {
            const double ratio=(1./result.temperatures[i-1]-1./result.temperatures[i])/(1./current[i-1]-1./current[i]);
            check(ratio>=std::exp(-.5)-1.e-8 && ratio<=std::exp(.5)+1.e-8,"inverse-temperature gap bound violated");
        }
        current=result.temperatures;
    }
    for (const auto& r : {std::vector<double>{-1.,.5},std::vector<double>{NAN,.5},std::vector<double>{1.,2.}}) {
        bool rejected=false;try{BATPT::UpdateSchedule(t,r,.5);}catch(const std::invalid_argument&){rejected=true;}
        check(rejected,"invalid rejection estimate accepted");
    }
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    try {
        BCLog::SetLogLevelScreen(BCLog::nothing); Model worker;
        if (rank) {
            double recv[3],value;
            for (;;) {
                MPI_Scatter(nullptr,3,MPI_DOUBLE,recv,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
                if(recv[0]==-1.)break;
                value=recv[0]==1.?worker.LogEval(std::vector<double>(recv+1,recv+3)):-INFINITY;
                MPI_Gather(&value,1,MPI_DOUBLE,nullptr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            }
        } else {
            schedule_tests();
            const bool joint=argc>1 && std::string(argv[1])=="joint";
            const std::string prefix=joint?"tune_joint_":"tune_factorized_";
            Model model; configure(model,joint);
            model.SetPTDiagnosticsOutput(prefix,73,5);
            model.SetProposalFunctionOutput(prefix+"proposal.txt");
            model.MarginalizeAll(BCIntegrate::kMargMetropolis);
            const auto status=model.GetPreRunStatus();
            check(status.temperature_adaptation && status.temperature_updates==3 && status.last_temperature_update==350 && status.temperatures_frozen,
                  "incorrect adaptation schedule/status");
            check(status.adaptation_sweeps==600 && status.settling_sweeps==31,"warmup budget changed");
            check(model.GetTemperatures().front()==1. && model.GetTemperatures().back()==16. &&
                  model.GetTemperatures()!=std::vector<double>({1.,1.01,2.,16.}),"ladder did not adapt");
            check(model.GetNChains()==8 && model.GetNColdChains()==2,"replica layout changed");
            check(std::abs(model.sum/model.draws-4./4.25)<.10,"cold mean incompatible with analytic posterior");
            const double variance=model.sum2/model.draws-std::pow(model.sum/model.draws,2);
            check(std::abs(variance-1./4.25)<.07,"cold variance incompatible with analytic posterior");
            check(model.GetStatistics().n_samples==10000,"warmup or hot samples contaminated production");
            // Each committed block uses one grid. No checkpoint may straddle an update.
            std::ifstream ledger(prefix+"pt_checkpoints.csv");std::string line;std::getline(ledger,line);
            while(std::getline(ledger,line)) {
                std::istringstream row(line);std::vector<std::string> cells;std::string cell;
                while(std::getline(row,cell,','))cells.push_back(cell);
                if(cells[0]!="prerun")continue;
                const unsigned first=std::stoul(cells[2]),last=std::stoul(cells[3]);
                for(unsigned update : {50u,150u,350u})check(!(first<=update && last>update),"diagnostic block straddles ladder change");
            }
            // A v2 cache restores its optimized intermediates, even with a different
            // configured initial grid.
            Model loaded;configure(loaded,joint);loaded.SetPTTemperatureAdaptation(false);
            loaded.SetParallelTempering({1.,3.,7.,16.},2,2);
            const std::string reload_prefix=prefix+"reload_";
            { std::ofstream stale(reload_prefix+"pt_temperature_adaptation.csv"); stale<<"stale history\n"; }
            loaded.SetPTDiagnosticsOutput(reload_prefix,73,0);
            loaded.LoadProposalFunction(prefix+"proposal.txt",7);
            loaded.SetNIterationsRun(19);loaded.MarginalizeAll(BCIntegrate::kMargMetropolis);
            check(loaded.GetTemperatures()==model.GetTemperatures() && loaded.kernel()==model.kernel(),"optimized cache did not restore matching grid and kernel");
            check(loaded.GetPreRunStatus().temperature_updates==0 && loaded.GetPreRunStatus().proposal_loaded,"reload repeated tuning");
            check(!std::ifstream(reload_prefix+"pt_temperature_adaptation.csv"),"reload retained stale temperature adaptation history");
            loaded.SaveProposalFunction(prefix+"restored.txt");
            // Validate invalid combinations before any MPI evaluation or initialization.
            for(unsigned problem=0;problem<8;++problem) {
                Model bad;configure(bad,joint);
                if(problem==0)bad.SetPreRunPolicy(BCEngineMCMC::kPreRunLegacy);
                if(problem==1)bad.SetParallelTempering({1.,16.},2,2);
                if(problem==2)bad.SetPreRunPolicy(BCEngineMCMC::kPreRunStaged,0);
                if(problem==3)bad.SetNIterationsPreRunMax(369);
                if(problem==4)bad.LoadProposalFunction(prefix+"proposal.txt",0);
                if(problem==5)bad.SetPTTemperatureAdaptation(true,1,3,.5);
                if(problem==6)bad.SetPTTemperatureAdaptation(true,1000000000,30,.5);
                if(problem==7)bad.SetNIterationsPreRunCheck(0);
                bool rejected=false;try{bad.MetropolisPreRun();}catch(const std::invalid_argument&){rejected=true;}
                check(rejected && bad.untouched(),"invalid PT adaptation reached initialization");
            }
            Model wrong;configure(wrong,joint);wrong.SetPTTemperatureAdaptation(false);
            wrong.SetParallelTempering({1.,2.,3.,32.},2,2);wrong.LoadProposalFunction(prefix+"proposal.txt",0);
            bool rejected=false;try{wrong.MetropolisPreRun();}catch(const std::runtime_error&){rejected=true;}
            check(rejected && wrong.untouched(),"optimized cache accepted different endpoints");
            bool late=false;try{model.SetPTTemperatureAdaptation(false);}catch(const std::logic_error&){late=true;}
            check(late,"late reconfiguration accepted");
            std::vector<double> stop(size*3,0.);double recv[3];for(int i=1;i<size;++i)stop[i*3]=-1.;
            if(size>1)MPI_Scatter(stop.data(),3,MPI_DOUBLE,recv,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
            std::cout<<"PASS: bounded barrier inversion, MPI temperature tuning, frozen kernels, analytic posterior, grid-aware cache, configuration validation\n";
        }
    } catch(const std::exception& e) {std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
    MPI_Finalize();
}
