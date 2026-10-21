// Streaming PT diagnostics. All counters have size fixed by slots and partitions,
// never by run length. CSV blocks are committed only after every block file closes.
#include "BCEngineMCMC.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <stdexcept>

namespace {
std::string quote(const std::string& value)
{
    std::string result = "\"";
    for (char c : value) { if (c == '"') result += '"'; result += c; }
    return result + '"';
}
std::ofstream csv(const std::string& path, bool append = true)
{
    std::ofstream out(path.c_str(), append ? std::ios::app : std::ios::trunc);
    if (!out) throw std::runtime_error("Cannot open PT diagnostics: " + path);
    out.imbue(std::locale::classic());
    out << std::setprecision(17);
    return out;
}
void closeCSV(std::ofstream& out, const std::string& path)
{
    out.flush();
    if (!out) throw std::runtime_error("Cannot flush PT diagnostics: " + path);
    out.close();
    if (!out) throw std::runtime_error("Cannot close PT diagnostics: " + path);
}
void commitCSV(std::ofstream& out, const std::string& path)
{
    closeCSV(out, path + ".tmp");
    if (std::rename((path + ".tmp").c_str(), path.c_str()) != 0)
        throw std::runtime_error("Cannot replace PT diagnostics: " + path);
}
void header(const std::string& path, const std::string& columns)
{
    auto out = csv(path, false);
    out << columns << '\n';
    closeCSV(out, path);
}
}

void BCEngineMCMC::SetPTDiagnosticsOutput(const std::string& prefix, unsigned checkpoint_sweeps,
                                         unsigned trace_sweeps)
{
    if (fMCMCPhase != kUnsetPhase || !fMCMCStates.empty())
        throw std::logic_error("Configure PT diagnostics before MCMC initialization");
    if (!prefix.empty() && (!GetParallelTempering() || !checkpoint_sweeps))
        throw std::invalid_argument("PT diagnostics require PT and a positive checkpoint interval");
    fPTDiagnostics = PTDiagnosticState();
    fPTDiagnostics.prefix = prefix;
    fPTDiagnostics.checkpoint_sweeps = checkpoint_sweeps;
    fPTDiagnostics.trace_sweeps = trace_sweeps;
}

void BCEngineMCMC::ClassifyPTDiagnosticState(const std::vector<double>& parameters,
                                            std::vector<unsigned>& labels) const
{
    if (fPTDiagnostics.partitions.empty()) { labels.clear(); return; }
    ClassifyPTDiagnosticRegions(parameters, labels);
    if (labels.size() != fPTDiagnostics.partitions.size())
        throw std::runtime_error("PT region classifier returned the wrong number of labels");
    for (unsigned p = 0; p < labels.size(); ++p)
        if (labels[p] >= fPTDiagnostics.partitions[p].labels.size())
            throw std::runtime_error("PT region classifier returned an invalid label index");
}

void BCEngineMCMC::BeginPTDiagnostics(const std::string& phase)
{
    auto& d = fPTDiagnostics;
    if (d.prefix.empty() || !GetParallelTempering()) return;
    if (d.active) EndPTDiagnostics();
    // A second sampling run on this engine replaces its previous output. Keep
    // prerun/warmup only when transitioning into their first production phase.
    if (d.initialized && (phase != "production" || d.phase == "production")) {
        const std::string prefix = d.prefix;
        const unsigned checkpoint = d.checkpoint_sweeps, trace = d.trace_sweeps;
        d = PTDiagnosticState();
        d.prefix = prefix; d.checkpoint_sweeps = checkpoint; d.trace_sweeps = trace;
    }
    if (!d.initialized) {
        d.partitions = DefinePTDiagnosticPartitions();
        std::set<std::string> names;
        for (const auto& p : d.partitions) {
            if (p.name.empty() || !names.insert(p.name).second || p.labels.empty())
                throw std::invalid_argument("PT partitions need unique names and nonempty label lists");
            std::set<std::string> labels;
            for (const auto& label : p.labels)
                if (label.empty() || !labels.insert(label).second)
                    throw std::invalid_argument("PT region labels must be nonempty and unique");
            const size_t r = p.labels.size();
            if (r > 1024 || d.nregions > std::numeric_limits<unsigned>::max() - r ||
                d.ntransitions > std::numeric_limits<unsigned>::max() - r*r)
                throw std::invalid_argument("PT partition schema is too large");
            d.region_offsets.push_back(d.nregions);
            d.transition_offsets.push_back(d.ntransitions);
            d.nregions += r;
            d.ntransitions += r*r;
        }
        header(d.prefix + "pt_local_moves.csv", "phase,block,first_sweep,last_sweep,ladder,slot,temperature,kernel,attempts,accepted,invalid_proposals,nonfinite_evaluations");
        header(d.prefix + "pt_region_occupancy.csv", "phase,block,first_sweep,last_sweep,ladder,slot,temperature,partition,region,samples,draws,weight");
        header(d.prefix + "pt_region_transitions.csv", "phase,block,first_sweep,last_sweep,ladder,slot,temperature,partition,from_region,to_region,cause,transitions");
        header(d.prefix + "pt_swaps_by_block.csv", "phase,block,first_sweep,last_sweep,ladder,temperature_low,temperature_high,attempts,accepted,acceptance_rate");
        header(d.prefix + "pt_roundtrip_regions.csv", "phase,block,first_sweep,last_sweep,ladder,replica,partition,departure_region,return_region,round_trips");
        header(d.prefix + "pt_checkpoints.csv", "phase,block,first_sweep,last_sweep");
        std::string columns = "phase,sweep,ladder,slot,temperature,replica,log_likelihood,log_prior";
        for (unsigned p = 0; p < d.partitions.size(); ++p) columns += ",region_" + std::to_string(p);
        header(d.prefix + "pt_states.csv", columns);
        auto schema = csv(d.prefix + "pt_partitions.csv", false);
        schema << "partition_index,partition,region_index,region\n";
        for (unsigned p = 0; p < d.partitions.size(); ++p)
            for (unsigned r = 0; r < d.partitions[p].labels.size(); ++r)
                schema << p << ',' << quote(d.partitions[p].name) << ',' << r << ',' << quote(d.partitions[p].labels[r]) << '\n';
        closeCSV(schema, d.prefix + "pt_partitions.csv");
        // Remove stale phase completion records when a prefix is reused for a new run.
        for (const char* old : {"prerun", "warmup", "production"})
            std::remove((d.prefix + "pt_progress_" + old + ".csv").c_str());
        d.initialized = true;
    }
    d.phase = phase;
    d.sweep = d.window_start = d.block = 0;
    d.slots.assign(GetNChains(), PTDiagnosticSlot());
    for (unsigned c = 0; c < GetNChains(); ++c) {
        auto& slot = d.slots[c];
        ClassifyPTDiagnosticState(fMCMCStates[c].parameters, slot.regions);
        slot.occupancy.assign(d.nregions, 0);
        slot.local_changes.assign(d.ntransitions, 0);
        slot.swap_changes.assign(d.ntransitions, 0);
    }
    d.swap_attempts.assign(fPT.ladder_attempts.size(), 0);
    d.swap_accepts.assign(fPT.ladder_accepts.size(), 0);
    d.trip_regions.assign(size_t(GetNChains()) * d.ntransitions, 0);
    d.departure_regions.assign(size_t(GetNChains()) * d.partitions.size(), 0);
    d.scratch.resize(d.partitions.size());
    for (unsigned c = 0; c < GetNChains(); ++c)
        for (unsigned p = 0; p < d.partitions.size(); ++p)
            d.departure_regions[size_t(fPT.replica_ids[c])*d.partitions.size()+p] = d.slots[c].regions[p];
    d.active = true;
    WritePTDiagnosticTrace();
    FlushPTDiagnostics(false); // Initial progress exists even before the first checkpoint.
}

void BCEngineMCMC::RecordPTDiagnosticLocalMove(unsigned c, bool multivariate, bool valid,
                                               bool finite, bool accepted)
{
    auto& d = fPTDiagnostics;
    if (!d.active) return;
    auto& s = d.slots[c];
    const unsigned kernel = multivariate ? 1 : 0;
    ++s.attempts[kernel];
    if (!valid) ++s.invalid[kernel];
    else if (!finite) ++s.nonfinite[kernel];
    if (!accepted) return;
    ++s.accepted[kernel];
    ClassifyPTDiagnosticState(fMCMCStates[c].parameters, d.scratch);
    for (unsigned p = 0; p < d.partitions.size(); ++p) {
        const unsigned from = s.regions[p], to = d.scratch[p], r = d.partitions[p].labels.size();
        if (from != to) ++s.local_changes[d.transition_offsets[p] + from*r + to];
        s.regions[p] = to;
    }
}

void BCEngineMCMC::RecordPTDiagnosticSwap(unsigned a, unsigned b, bool accepted)
{
    auto& d = fPTDiagnostics;
    if (!d.active) return;
    ++d.swap_attempts[a]; // temperature-major low-slot index equals edge*nladders+ladder
    if (!accepted) return;
    ++d.swap_accepts[a];
    const unsigned np = d.partitions.size();
    if (a < GetNColdChains()) // Record the actual cold departure, after any local moves.
        for (unsigned p = 0; p < np; ++p)
            d.departure_regions[size_t(fPT.replica_ids[a])*np+p] = d.slots[a].regions[p];
    for (unsigned p = 0; p < np; ++p) {
        const unsigned ra = d.slots[a].regions[p], rb = d.slots[b].regions[p], r = d.partitions[p].labels.size();
        if (ra != rb) {
            ++d.slots[a].swap_changes[d.transition_offsets[p]+ra*r+rb];
            ++d.slots[b].swap_changes[d.transition_offsets[p]+rb*r+ra];
        }
        std::swap(d.slots[a].regions[p], d.slots[b].regions[p]);
    }
}

void BCEngineMCMC::RecordPTDiagnosticRoundTrip(unsigned replica, unsigned cold_slot)
{
    auto& d = fPTDiagnostics;
    if (!d.active) return;
    for (unsigned p = 0; p < d.partitions.size(); ++p) {
        const unsigned from = d.departure_regions[size_t(replica)*d.partitions.size()+p];
        const unsigned to = d.slots[cold_slot].regions[p], r = d.partitions[p].labels.size();
        ++d.trip_regions[size_t(replica)*d.ntransitions+d.transition_offsets[p]+from*r+to];
    }
}

void BCEngineMCMC::WritePTDiagnosticTrace()
{
    const auto& d = fPTDiagnostics;
    if (!d.active || !d.trace_sweeps) return;
    auto out = csv(d.prefix + "pt_states.csv");
    for (unsigned c = 0; c < GetNChains(); ++c) {
        out << d.phase << ',' << d.sweep << ',' << c%GetNColdChains() << ',' << c << ','
            << GetChainTemperature(c) << ',' << fPT.replica_ids[c] << ','
            << fMCMCStates[c].log_likelihood << ',' << fMCMCStates[c].log_prior;
        for (unsigned r : d.slots[c].regions) out << ',' << r;
        out << '\n';
    }
    closeCSV(out, d.prefix + "pt_states.csv");
}

void BCEngineMCMC::FinishPTDiagnosticsSweep()
{
    auto& d = fPTDiagnostics;
    if (!d.active) return;
    d.sweep = fMCMCCurrentIteration;
    for (auto& slot : d.slots)
        for (unsigned p = 0; p < d.partitions.size(); ++p)
            ++slot.occupancy[d.region_offsets[p] + slot.regions[p]];
    if (d.trace_sweeps && d.sweep % d.trace_sweeps == 0) WritePTDiagnosticTrace();
    if (d.sweep % d.checkpoint_sweeps == 0) FlushPTDiagnostics(false);
}

void BCEngineMCMC::EndPTDiagnostics()
{
    auto& d = fPTDiagnostics;
    if (!d.active) return;
    if (d.trace_sweeps && d.sweep && d.sweep % d.trace_sweeps != 0) WritePTDiagnosticTrace();
    FlushPTDiagnostics(true);
    d.active = false;
}

void BCEngineMCMC::FlushPTDiagnostics(bool complete)
{
    auto& d = fPTDiagnostics;
    if (!d.active) return;
    const unsigned nl = GetNColdChains();
    if (d.sweep > d.window_start) {
        auto moves = csv(d.prefix + "pt_local_moves.csv");
        auto occupancy = csv(d.prefix + "pt_region_occupancy.csv");
        auto changes = csv(d.prefix + "pt_region_transitions.csv");
        auto swaps = csv(d.prefix + "pt_swaps_by_block.csv");
        auto trips = csv(d.prefix + "pt_roundtrip_regions.csv");
        auto prefix = [&](std::ofstream& out) -> std::ofstream& {
            out << d.phase << ',' << d.block << ',' << d.window_start+1 << ',' << d.sweep << ',';
            return out;
        };
        for (unsigned c = 0; c < GetNChains(); ++c) {
            const auto& s = d.slots[c];
            for (unsigned k = 0; k < 2; ++k) if (s.attempts[k])
                prefix(moves) << c%nl << ',' << c << ',' << GetChainTemperature(c) << ','
                    << (k ? "multivariate" : "factorized") << ',' << s.attempts[k] << ',' << s.accepted[k]
                    << ',' << s.invalid[k] << ',' << s.nonfinite[k] << '\n';
            for (unsigned p = 0; p < d.partitions.size(); ++p) {
                const auto& part = d.partitions[p];
                const unsigned r = part.labels.size();
                for (unsigned from = 0; from < r; ++from) {
                    const auto count = s.occupancy[d.region_offsets[p]+from];
                    prefix(occupancy) << c%nl << ',' << c << ',' << GetChainTemperature(c) << ',' << quote(part.name)
                        << ',' << quote(part.labels[from]) << ',' << count << ',' << d.sweep-d.window_start
                        << ',' << double(count)/(d.sweep-d.window_start) << '\n';
                    for (unsigned to = 0; to < r; ++to) {
                        const unsigned index = d.transition_offsets[p]+from*r+to;
                        const unsigned long long counts[] = {s.local_changes[index], s.swap_changes[index]};
                        for (unsigned cause = 0; cause < 2; ++cause) if (counts[cause])
                            prefix(changes) << c%nl << ',' << c << ',' << GetChainTemperature(c) << ',' << quote(part.name)
                                << ',' << quote(part.labels[from]) << ',' << quote(part.labels[to]) << ','
                                << (cause ? "swap" : "local") << ',' << counts[cause] << '\n';
                        const auto ntrip = d.trip_regions[size_t(c)*d.ntransitions+index];
                        if (ntrip) prefix(trips) << c%nl << ',' << c << ',' << quote(part.name) << ','
                            << quote(part.labels[from]) << ',' << quote(part.labels[to]) << ',' << ntrip << '\n';
                    }
                }
            }
        }
        for (unsigned t = 0; t+1 < fPT.temperatures.size(); ++t)
            for (unsigned l = 0; l < nl; ++l) {
                const unsigned k = t*nl+l;
                prefix(swaps) << l << ',' << fPT.temperatures[t] << ',' << fPT.temperatures[t+1] << ','
                    << d.swap_attempts[k] << ',' << d.swap_accepts[k] << ',';
                if (d.swap_attempts[k]) swaps << double(d.swap_accepts[k])/d.swap_attempts[k];
                else swaps << "nan";
                swaps << '\n';
            }
        closeCSV(moves, d.prefix + "pt_local_moves.csv");
        closeCSV(occupancy, d.prefix + "pt_region_occupancy.csv");
        closeCSV(changes, d.prefix + "pt_region_transitions.csv");
        closeCSV(swaps, d.prefix + "pt_swaps_by_block.csv");
        closeCSV(trips, d.prefix + "pt_roundtrip_regions.csv");
        // The ledger is the commit record. Readers ignore rows beyond this boundary.
        auto ledger = csv(d.prefix + "pt_checkpoints.csv");
        ledger << d.phase << ',' << d.block << ',' << d.window_start+1 << ',' << d.sweep << '\n';
        closeCSV(ledger, d.prefix + "pt_checkpoints.csv");
        for (auto& slot : d.slots) {
            slot.ClearMoves();
            std::fill(slot.occupancy.begin(), slot.occupancy.end(), 0);
            std::fill(slot.local_changes.begin(), slot.local_changes.end(), 0);
            std::fill(slot.swap_changes.begin(), slot.swap_changes.end(), 0);
        }
        std::fill(d.swap_attempts.begin(), d.swap_attempts.end(), 0);
        std::fill(d.swap_accepts.begin(), d.swap_accepts.end(), 0);
        std::fill(d.trip_regions.begin(), d.trip_regions.end(), 0);
        d.window_start = d.sweep;
        ++d.block;
    }
    if (d.phase == "production") {
        WriteParallelTemperingDiagnostics(d.prefix + "pt_swaps.csv");
        WriteParallelTemperingTransport(d.prefix);
    }
    const std::string path = d.prefix + "pt_progress_" + d.phase + ".csv";
    auto progress = csv(path + ".tmp", false);
    progress << "phase,last_sweep,completed,checkpoint_sweeps,trace_sweeps\n" << d.phase << ',' << d.sweep << ','
             << (complete ? 1 : 0) << ',' << d.checkpoint_sweeps << ',' << d.trace_sweeps << '\n';
    commitCSV(progress, path);
}

void BCEngineMCMC::WriteParallelTemperingDiagnostics(const std::string& filename) const
{
    auto out = csv(filename + ".tmp", false);
    if (!out) throw std::runtime_error("Cannot write PT diagnostics: " + filename);
    out << std::setprecision(17) << "temperature_low,temperature_high,attempts,accepted,acceptance_rate\n";
    for (unsigned t = 0; t < fPT.attempts.size(); ++t)
        out << fPT.temperatures[t] << ',' << fPT.temperatures[t+1] << ','
            << fPT.attempts[t] << ',' << fPT.accepts[t] << ','
            << (fPT.attempts[t] ? double(fPT.accepts[t]) / fPT.attempts[t] : 0.) << '\n';
    commitCSV(out, filename);
}

void BCEngineMCMC::WriteParallelTemperingTransport(const std::string& prefix) const
{
    if (!GetParallelTempering())
        throw std::logic_error("Replica transport diagnostics require parallel tempering");
    auto trips = csv(prefix + "pt_roundtrips.csv.tmp", false);
    auto swaps = csv(prefix + "pt_swaps_by_ladder.csv.tmp", false);
    if (!trips || !swaps) throw std::runtime_error("Cannot write PT transport diagnostics: " + prefix);
    trips << std::setprecision(17);
    swaps << std::setprecision(17);
    trips << "replica,ladder,round_trips,total_sweeps,mean_sweeps,min_sweeps,max_sweeps,journey_stage\n";
    const unsigned ladders = GetNColdChains();
    for (unsigned id = 0; id < fPT.round_trips.size(); ++id) {
        trips << id << ',' << id % ladders << ',' << fPT.round_trips[id] << ',' << fPT.trip_sweeps[id] << ',';
        if (fPT.round_trips[id]) trips << double(fPT.trip_sweeps[id]) / fPT.round_trips[id] << ','
            << fPT.min_trip_sweeps[id] << ',' << fPT.max_trip_sweeps[id];
        else trips << "nan,nan,nan";
        trips << ',' << fPT.journey_stage[id] << '\n';
    }
    swaps << "ladder,temperature_low,temperature_high,attempts,accepted,acceptance_rate\n";
    for (unsigned t = 0; t < fPT.attempts.size(); ++t)
        for (unsigned l = 0; l < ladders; ++l) {
            const unsigned k = t*ladders+l;
            swaps << l << ',' << fPT.temperatures[t] << ',' << fPT.temperatures[t+1] << ','
                << fPT.ladder_attempts[k] << ',' << fPT.ladder_accepts[k] << ',';
            if (fPT.ladder_attempts[k]) swaps << double(fPT.ladder_accepts[k]) / fPT.ladder_attempts[k];
            else swaps << "nan";
            swaps << '\n';
        }
    commitCSV(trips, prefix + "pt_roundtrips.csv");
    commitCSV(swaps, prefix + "pt_swaps_by_ladder.csv");
}

