// Numerical guards and end-to-end distribution checks for smooth scale tuning.
#include <BAT/BCModel.h>
#include <BAT/BCLog.h>
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
static void check(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
static bool near(double a,double b) {return std::abs(a-b)<1.e-12*std::max(1.,std::abs(b));}
class Probe : public BCModel {
public:
    Probe() : BCModel("smooth_probe") {
        AddParameter("x",-8.,8.);AddParameter("fixed",-1.,1.);AddParameter("y",-8.,8.);
        SetPriorConstantAll();GetParameter(1).Fix(.4);
    }
    double LogLikelihood(const std::vector<double>& p) override {return -.5*(p[0]*p[0]+p[2]*p[2]);}
    void prepare(double gain=2.,double decay=.6) {
        SetScaleAdaptation(kScaleSmooth,.25,gain,decay);SetNChains(2);MCMCInitialize();
        fMCMCPhase=kPreRun;fWarmup.status.stage="adapting";
        fMCMCProposalFunctionScaleFactor.assign(2,std::vector<double>{1.,0.,1.});
    }
    double step(double rate,bool joint=false) {
        for(auto& s:fMCMCStatistics) {s.n_samples_efficiency=100;s.efficiency.assign(3,rate);}
        UpdateSmoothProposalScales(joint);return fMCMCProposalFunctionScaleFactor[0][0];
    }
    double fixedScale() const {return fMCMCProposalFunctionScaleFactor[0][1];}
    void freeze() {fWarmup.fixed_kernel=true;}
};
class Gaussian : public BCModel {
public:
    unsigned draws=0;double sum=0.,squares=0.;
    Gaussian() : BCModel("smooth_gaussian") {
        AddParameter("fixed",-1.,1.);AddParameter("x",-8.,8.);AddParameter("y",-8.,8.);
        SetPriorConstantAll();GetParameter(0).Fix(.4);
    }
    double LogLikelihood(const std::vector<double>& p) override {return -.5*(p[1]*p[1]+p[2]*p[2]);}
    void MCMCUserIterationInterface() override {
        for(unsigned c=0;c<GetNColdChains();++c) {const double x=fMCMCStates[c].parameters[1];sum+=x;squares+=x*x;++draws;}
    }
};
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    try {
        BCLog::SetLogLevelScreen(BCLog::nothing);
        Probe up;up.prepare();check(near(up.step(1.),std::exp(.5)),"upper log-step cap or positive update wrong");
        check(up.fixedScale()==0.,"fixed coordinate scale modified");
        Probe down;down.prepare();check(near(down.step(0.),std::exp(-.5)),"negative update wrong");
        Probe decay;decay.prepare();const double a=decay.step(.3),b=decay.step(.3);
        check(near(std::log(a),.1) && near(std::log(b/a),.1/std::pow(2.,.6)),"learning rate does not diminish as specified");
        Probe joint;joint.prepare();check(near(joint.step(.3,true),std::exp(.2)),"joint multiplier was treated as standard deviation instead of variance");
        check(joint.GetPreRunStatus().joint_scale_updates==1 && joint.GetPreRunStatus().factorized_scale_updates==0,"kernel gain counters conflated");
        Probe center;center.prepare();check(near(center.step(.25),1.),"scale changed at target acceptance");
        Probe bounds;bounds.prepare();bounds.SetScaleFactorLowerLimit(.9);bounds.SetScaleFactorUpperLimit(1.1);
        check(near(bounds.step(1.),1.1) && near(bounds.step(0.),.9),"scale limits not respected");
        center.freeze();bool rejected=false;try{center.step(.3);}catch(const std::logic_error&){rejected=true;}
        check(rejected,"adaptation occurred after freezing");
        for(unsigned bad=0;bad<6;++bad) {
            Probe p;bool failed=false;
            try {p.SetScaleAdaptation(BCEngineMCMC::kScaleSmooth,bad==0?0.:bad==1?1.:.25,
                    bad==2?0.:bad==3?INFINITY:2.,bad==4?.5:bad==5?1.1:.6);}catch(const std::invalid_argument&){failed=true;}
            check(failed,"invalid tuning hyperparameter accepted");
        }
        for(unsigned bad=0;bad<3;++bad) {
            Probe p;p.SetScaleAdaptation(BCEngineMCMC::kScaleSmooth);
            if(bad==0)p.SetScaleFactorLowerLimit(-1.);
            if(bad==1)p.SetScaleFactorUpperLimit(0.);
            if(bad==2){p.SetScaleFactorLowerLimit(2.);p.SetScaleFactorUpperLimit(1.);}
            bool failed=false;try{p.MetropolisPreRun();}catch(const std::invalid_argument&){failed=true;}
            check(failed,"invalid scale limits accepted");
        }
        for(bool multivariate:{false,true})for(bool pt:{false,true}) {
            Gaussian m;m.SetNChains(2);if(pt)m.SetParallelTempering({1.,2.,4.},2,2);
            m.SetProposeMultivariate(multivariate);m.SetProposalFunctionDof(0.);
            m.SetScaleAdaptation(BCEngineMCMC::kScaleSmooth);m.SetPreRunPolicy(BCEngineMCMC::kPreRunStaged,200);
            m.SetNIterationsPreRunMax(2000);m.SetNIterationsPreRunCheck(100);m.SetNIterationsPreRunFactorized(500);
            m.SetNIterationsRun(12000);m.SetInitialPositionScheme(BCEngineMCMC::kInitRandomPrior);m.SetRandomSeed(89437);
            m.MarginalizeAll(BCIntegrate::kMargMetropolis);
            const double mean=m.sum/m.draws,var=m.squares/m.draws-mean*mean;
            std::cout<<"joint="<<multivariate<<" PT="<<pt<<" cold mean="<<mean<<" variance="<<var<<'\n';
            check(std::abs(mean)<.1 && std::abs(var-1.)<.15,"smooth tuning changed the Gaussian posterior");
            check(m.draws==24000 && m.GetStatistics().n_samples==24000,"hot or discarded samples contaminated posterior");
            check(m.GetPreRunStatus().factorized_scale_updates>0 && (!multivariate || m.GetPreRunStatus().joint_scale_updates>0),"gain schedules missing");
        }
        std::cout<<"PASS: signed bounded updates, diminishing gain, variance units, fixed coordinates, invalid settings, Gaussian posterior\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
    MPI_Finalize();
}
