// Event-level tests: distinguish local moves from swaps, count brief crossings,
// and keep a committed checkpoint usable without a normal end-of-run call.
#include <BAT/BCModel.h>
#include <BAT/BCLog.h>
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

static void check(bool condition, const char* reason) { if (!condition) throw std::runtime_error(reason); }
using Row = std::map<std::string, std::string>;
static std::vector<Row> readCSV(const std::string& path) {
    std::ifstream input(path); check(bool(input), "CSV missing");
    std::string line, item; std::vector<std::string> names; std::vector<Row> rows;
    std::getline(input,line); std::istringstream h(line);
    while(std::getline(h,item,',')) names.push_back(item);
    while(std::getline(input,line)) {
        Row row; std::istringstream fields(line); unsigned i=0;
        while(std::getline(fields,item,',')) {
            if(item.size()>1 && item.front()=='"' && item.back()=='"') item=item.substr(1,item.size()-2);
            row[names.at(i++)]=item;
        }
        rows.push_back(row);
    }
    return rows;
}
class FlatDiagnostic : public BCModel {
public:
    FlatDiagnostic() : BCModel("stream_test") { AddParameter("x",-1.,1.); SetPriorConstantAll(); }
    double LogLikelihood(const std::vector<double>&) override { return 0.; }
    std::vector<PTDiagnosticPartition> DefinePTDiagnosticPartitions() const override {
        return {{"sign", {"negative", "positive"}}};
    }
    void ClassifyPTDiagnosticRegions(const std::vector<double>& p, std::vector<unsigned>& labels) const override {
        labels.assign(1,p.at(0)<0.?0:1);
    }
    void start(const std::string& prefix, unsigned trace=2) {
        SetParallelTempering({1.,2.},1,1); SetPTDiagnosticsOutput(prefix,2,trace);
        SetInitialPositionScheme(kInitCenter); MCMCInitialize();
        fMCMCStates[0].parameters[0]=-.5; fMCMCStates[1].parameters[0]=.5;
        fMCMCPhase=kPreRun; fMCMCCurrentIteration=0; BeginPTDiagnostics("prerun");
    }
    void local(unsigned slot, double x, bool accepted=true, bool valid=true, bool finite=true) {
        if(accepted) fMCMCStates[slot].parameters[0]=x;
        RecordPTDiagnosticLocalMove(slot,false,valid,finite,accepted);
    }
    void sweep(bool exchange=false) { ++fMCMCCurrentIteration; if(exchange) ExchangeReplicas(); FinishPTDiagnosticsSweep(); }
    void production() {
        EndPTDiagnostics(); fMCMCPhase=kMainRun; fMCMCCurrentIteration=0;
        ResetReplicaExchangeStatistics(); BeginPTDiagnostics("production");
    }
    void end() { EndPTDiagnostics(); }
    void rejectNextSwap() { fMCMCStates[0].log_likelihood=0.; fMCMCStates[1].log_likelihood=-1.e300; }
    size_t cells() const {
        size_t n=fPTDiagnostics.trip_regions.size()+fPTDiagnostics.departure_regions.size();
        for(const auto& s:fPTDiagnostics.slots) n+=s.regions.size()+s.occupancy.size()+s.local_changes.size()+s.swap_changes.size();
        return n;
    }
};
class SampledDiagnostic : public FlatDiagnostic {
public:
    std::vector<double> history;
    void MCMCUserIterationInterface() override {
        for (unsigned c=0;c<GetNChains();++c) {
            history.push_back(fMCMCStates[c].parameters[0]);
            history.push_back(GetReplicaIDs()[c]);
        }
    }
    void run(bool diagnostics, bool joint) {
        SetParallelTempering({1.,2.,4.},2,2);
        if(diagnostics) SetPTDiagnosticsOutput(joint ? "neutral_joint_" : "neutral_factorized_",13,7);
        SetProposeMultivariate(joint); SetProposalFunctionDof(0.);
        SetNIterationsPreRunMin(100); SetNIterationsPreRunMax(200);
        SetNIterationsPreRunFactorized(50); SetNIterationsPreRunCheck(50);
        SetNIterationsRun(300); SetInitialPositionScheme(kInitRandomPrior);
        SetRandomSeed(62733);
        MarginalizeAll(BCIntegrate::kMargMetropolis);
    }
};
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    try {
        BCLog::SetLogLevelScreen(BCLog::nothing);
        for(bool joint : {false,true}) {
            SampledDiagnostic baseline, observed;
            baseline.run(false,joint); observed.run(true,joint);
            check(baseline.history==observed.history,"diagnostics changed the sampled path or replica IDs");
        }
        FlatDiagnostic m; m.start("stream_");
        m.local(0,.5); m.sweep(); m.local(0,-.5); m.sweep(); m.sweep();
        m.production();
        check(readCSV("stream_pt_progress_prerun.csv")[0]["completed"]=="1", "prerun not finalized");
        m.local(0,-.5,false,false); m.local(1,.7); m.sweep(true);
        // Three local crossings occur between trace samples. They must all be counted.
        m.local(1,.8); m.local(1,-.7); m.local(1,.8); m.sweep(true);
        auto regions=readCSV("stream_pt_region_transitions.csv");
        unsigned local_hot=0,swap_cold=0;
        for(const auto& r:regions) if(r.at("phase")=="production") {
            if(r.at("slot")=="1" && r.at("cause")=="local") local_hot+=std::stoul(r.at("transitions"));
            if(r.at("slot")=="0" && r.at("cause")=="swap") swap_cold+=std::stoul(r.at("transitions"));
        }
        check(local_hot==3 && swap_cold==1,"crossings were aliased or local and swap events were conflated");
        auto trips=readCSV("stream_pt_roundtrip_regions.csv");
        check(trips.size()==1 && trips[0]["replica"]=="0" && trips[0]["departure_region"]=="negative" &&
              trips[0]["return_region"]=="positive" && trips[0]["round_trips"]=="1", "mode-changing trip not counted");
        check(m.GetReplicaRoundTrips()[0]==1 && m.GetReplicaRoundTrips()[1]==0,"hot-origin replica counted too early");
        auto moves=readCSV("stream_pt_local_moves.csv");
        for(const auto& r:moves) if(r.at("phase")=="production" && r.at("slot")=="0")
            check(r.at("attempts")=="1" && r.at("accepted")=="0" && r.at("invalid_proposals")=="1", "invalid proposal miscounted");
        auto occupancy=readCSV("stream_pt_region_occupancy.csv");
        for(const auto& r:occupancy) if(r.at("phase")=="production" && r.at("slot")=="0" && r.at("region")=="positive")
            check(r.at("samples")=="2" && r.at("draws")=="2", "occupancy not sampled after exchanges");
        const auto cells=m.cells();
        for(unsigned i=0;i<10000;++i) m.local(0,i%2?.5:-.5);
        check(m.cells()==cells,"diagnostic storage grew with event count");
        // No end() here simulates interruption after one additional, uncommitted sweep.
        m.sweep();
        auto ledger=readCSV("stream_pt_checkpoints.csv");
        check(ledger.back()["phase"]=="production" && ledger.back()["last_sweep"]=="2", "checkpoint includes uncommitted data");
        auto progress=readCSV("stream_pt_progress_production.csv");
        check(progress[0]["last_sweep"]=="2" && progress[0]["completed"]=="0", "interruption not represented");
        check(!readCSV("stream_pt_swaps_by_ladder.csv").empty(),"periodic legacy transport snapshot absent");
        m.end();
        check(readCSV("stream_pt_progress_production.csv")[0]["completed"]=="1", "completion not recorded");
        FlatDiagnostic silent; silent.start("silent_",0); silent.sweep(); silent.end();
        check(readCSV("silent_pt_states.csv").empty(),"trace disabling ignored");
        FlatDiagnostic failedSwap; failedSwap.start("rejected_",0);
        failedSwap.rejectNextSwap(); failedSwap.sweep(true); failedSwap.end();
        const auto rejectedSwaps=readCSV("rejected_pt_swaps_by_block.csv");
        check(rejectedSwaps.size()==1 && rejectedSwaps[0].at("attempts")=="1" && rejectedSwaps[0].at("accepted")=="0",
              "rejected swap miscounted");
        check(readCSV("rejected_pt_region_transitions.csv").empty(),"rejected swap changed region counters");
        FlatDiagnostic reused; reused.start("reused_",0); reused.sweep(); reused.end();
        reused.production(); reused.sweep(); reused.end();
        reused.production(); reused.sweep(); reused.end();
        const auto reusedLedger=readCSV("reused_pt_checkpoints.csv");
        check(reusedLedger.size()==1 && reusedLedger[0].at("phase")=="production",
              "repeated engine run appended conflicting phase blocks");
        FlatDiagnostic invalid; bool rejected=false;
        try { invalid.SetPTDiagnosticsOutput("invalid_",0); } catch(const std::invalid_argument&) { rejected=true; }
        check(rejected,"invalid diagnostic configuration accepted");
        std::cout<<"PASS: local/swap attribution, exact crossings, changed-region trips, independent phases, interruption, bounded counters, trace disabling\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<std::endl; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Finalize();
}
