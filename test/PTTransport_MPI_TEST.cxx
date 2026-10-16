// Deterministic transport tests: flat targets make every exchange accepted.
#include <BAT/BCModel.h>
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

static void check(bool ok, const char* why) { if(!ok) throw std::runtime_error(why); }
class Flat : public BCModel {
public:
    Flat() : BCModel("transport_test") { AddParameter("x",-1.,1.); SetPriorConstantAll(); }
    double LogLikelihood(const std::vector<double>&) override { return 0.; }
    void start(const std::vector<double>& temperatures, unsigned ladders) {
        SetParallelTempering(temperatures,ladders,10);
        SetInitialPositionScheme(kInitCenter);
        MCMCInitialize(); fMCMCPhase=kMainRun; ResetReplicaExchangeStatistics();
    }
    void step(unsigned sweep) { fMCMCCurrentIteration=sweep; ExchangeReplicas(); }
    void reset() { ResetReplicaExchangeStatistics(); }
    void reject() { fMCMCStates[0].log_likelihood=100.; fMCMCStates[1].log_likelihood=-1e300; }
    unsigned long long duration(unsigned id) const { return fPT.trip_sweeps.at(id); }
    unsigned long long ladderAttempts(unsigned edge,unsigned ladder) const {
        return fPT.ladder_attempts.at(edge*GetNColdChains()+ladder);
    }
};
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    try {
        Flat two; two.start({1.,4.},1);
        check(two.GetReplicaRoundTrips()==std::vector<unsigned long long>({0,0}),"bad initial counts");
        two.step(10);
        check(two.GetReplicaIDs()==std::vector<unsigned>({1,0}),"IDs not exchanged");
        check(two.GetReplicaRoundTrips()==std::vector<unsigned long long>({0,0}),"hot-origin arrival counted as a round trip");
        two.step(20);
        check(two.GetReplicaRoundTrips()==std::vector<unsigned long long>({1,0}),"cold-origin round trip missing");
        check(two.duration(0)==20,"incorrect round-trip duration");
        two.step(30);
        check(two.GetReplicaRoundTrips()==std::vector<unsigned long long>({1,1}),"hot-origin trip needs a prior cold visit");
        two.WriteParallelTemperingTransport("transport_test_");
        std::ifstream csv("transport_test_pt_roundtrips.csv"); std::string header;
        std::getline(csv,header); check(header.find("mean_sweeps")!=std::string::npos,"missing transport output");
        two.reset();
        check(two.GetReplicaRoundTrips()==std::vector<unsigned long long>({0,0}),"warm-up trips leaked into main run");
        check(two.GetReplicaIDs()==std::vector<unsigned>({0,1}),"IDs not reset at production boundary");
        two.reject(); two.step(10);
        check(two.GetReplicaIDs()==std::vector<unsigned>({0,1}),"rejected exchange moved IDs");
        check(two.GetSwapAccepts()[0]==0 && two.GetSwapAttempts()[0]==1,"bad rejection counters");
        Flat three; three.start({1.,2.,4.},2);
        for(unsigned t=10;t<=50;t+=10) three.step(t);
        check(three.GetReplicaRoundTrips()[0]==1 && three.GetReplicaRoundTrips()[1]==1,"three-temperature trip missing");
        check(three.duration(0)==50,"alternating pairs duration incorrect");
        for(unsigned slot=0;slot<6;++slot) check(three.GetReplicaIDs()[slot]%2==slot%2,"configuration crossed independent ladders");
        for(unsigned edge=0;edge<2;++edge)
            check(three.ladderAttempts(edge,0)+three.ladderAttempts(edge,1)==three.GetSwapAttempts()[edge],"per-ladder counters do not sum to aggregate");
        std::cout<<"PASS transport: physical IDs, cold-hot-cold cycles, hot-origin exclusion, durations, reset, rejected swaps, independent ladders"<<std::endl;
    } catch(const std::exception& e) { std::cerr<<e.what()<<std::endl; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Finalize();
}
