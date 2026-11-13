// Prerun-only temperature schedule adaptation; no extra likelihood calls or RNG.
#include "BCEngineMCMC.h"
#include "PTTemperatureSchedule.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <stdexcept>

void BCEngineMCMC::SetPTTemperatureAdaptation(bool enabled, unsigned initial_sweeps,
                                             unsigned rounds, double damping)
{
    if (fMCMCPhase != kUnsetPhase || !fMCMCStates.empty())
        throw std::logic_error("Configure PT temperature adaptation before MCMC initialization");
    if (!initial_sweeps || !rounds || rounds > 30 || !std::isfinite(damping) || damping <= 0. || damping > 1.)
        throw std::invalid_argument("PT temperature adaptation needs positive batch/round counts (rounds <=30) and damping in (0,1]");
    fPTTuning.enabled = enabled; fPTTuning.initial_sweeps = initial_sweeps;
    fPTTuning.rounds = rounds; fPTTuning.damping = damping;
    ResetPreRunStatus();
}

void BCEngineMCMC::SetPTTemperatureAdaptationOutput(const std::string& filename)
{
    if (fMCMCPhase != kUnsetPhase || !fMCMCStates.empty())
        throw std::logic_error("Configure PT temperature history before MCMC initialization");
    fPTTuning.output = filename;
}

void BCEngineMCMC::ValidatePTTemperatureAdaptation() const
{
    if (!fPTTuning.enabled) return;
    if (!GetParallelTempering() || fPT.temperatures.size() < 3 || fWarmup.policy != kPreRunStaged ||
        !fWarmup.settling_sweeps || !fProposalCache.input.empty())
        throw std::invalid_argument("PT temperature adaptation requires >=3 temperatures, staged prerun, positive settling and no loaded proposal");
    const unsigned long long total = static_cast<unsigned long long>(fPTTuning.initial_sweeps)*
        ((1ULL << fPTTuning.rounds)-1);
    const unsigned long long last = std::max<unsigned long long>(total,
        fMCMCProposeMultivariate ? fMCMCNIterationsPreRunFactorized : 0);
    if (!fMCMCNIterationsPreRunCheck || last + fMCMCNIterationsPreRunCheck > fMCMCNIterationsPreRunMax)
        throw std::invalid_argument("PT temperature adaptation batches must leave at least one full proposal-check window on the final ladder (and joint kernel) within the prerun budget");
    if (fPTTuning.initial_sweeps < fPT.swap_interval)
        throw std::invalid_argument("PT temperature adaptation initial batch must include at least one exchange round");
    BATPT::UpdateSchedule(fPT.temperatures, std::vector<double>(fPT.temperatures.size()-1, .5), fPTTuning.damping);
}

void BCEngineMCMC::InitializePTTemperatureAdaptation()
{
    if (!fPTTuning.enabled) return;
    fPTTuning.next_sweep = fPTTuning.batch_sweeps = fPTTuning.initial_sweeps;
    fPTTuning.rejection_sum.assign(fPT.temperatures.size()-1, 0.);
    fPTTuning.observations.assign(fPT.temperatures.size()-1, 0);
    if (fPTTuning.output.empty() && !fPTDiagnostics.prefix.empty())
        fPTTuning.output = fPTDiagnostics.prefix + "pt_temperature_adaptation.csv";
    WritePTTemperatureUpdate(fPT.temperatures, {}, 0., true);
    BCLog::OutSummary(Form(" --> PT temperature adaptation: %u doubling batches from %u sweeps, damping=%g; fixed endpoints and replica count.",
        fPTTuning.rounds, fPTTuning.initial_sweeps, fPTTuning.damping));
}

void BCEngineMCMC::RecordPTTemperatureRejections()
{
    if (!fPTTuning.enabled || fMCMCPhase != kPreRun || fWarmup.fixed_kernel ||
        fWarmup.status.stage != "adapting" || fWarmup.status.temperatures_frozen) return;
    const unsigned ladders = GetNColdChains();
    // Measure ALL neighboring pairs before performing any disjoint swaps, as in
    // Syed et al., Algorithm 1. Expected rejection avoids Bernoulli noise.
    for (unsigned t = 0; t+1 < fPT.temperatures.size(); ++t)
        for (unsigned l = 0; l < ladders; ++l) {
            const double delta = (1./fPT.temperatures[t]-1./fPT.temperatures[t+1]) *
                (fMCMCStates[(t+1)*ladders+l].log_likelihood-fMCMCStates[t*ladders+l].log_likelihood);
            if (std::isnan(delta)) throw std::runtime_error("Nonfinite PT temperature adaptation swap ratio");
            fPTTuning.rejection_sum[t] += -std::expm1(std::min(0., delta));
            ++fPTTuning.observations[t];
        }
}

bool BCEngineMCMC::UpdatePTTemperatures()
{
    auto& status = fWarmup.status;
    if (!fPTTuning.enabled || status.temperatures_frozen ||
        unsigned(fMCMCCurrentIteration) < fPTTuning.next_sweep) return false;
    if (fMCMCPhase != kPreRun || fWarmup.fixed_kernel || status.stage != "adapting")
        throw std::logic_error("PT temperature updates are allowed only during adaptive prerun");
    std::vector<double> rejection(fPTTuning.rejection_sum.size());
    for (unsigned t = 0; t < rejection.size(); ++t) {
        if (!fPTTuning.observations[t]) throw std::runtime_error("PT temperature batch has no rejection observations");
        rejection[t] = fPTTuning.rejection_sum[t]/fPTTuning.observations[t];
    }
    const auto old = fPT.temperatures;
    const auto update = BATPT::UpdateSchedule(old, rejection, fPTTuning.damping);
    // A committed diagnostic block must never straddle a temperature change.
    FlushPTDiagnostics(false);
    fPT.temperatures = update.temperatures;
    fProposalCache.optimized_ladder = true;
    ++status.temperature_updates; status.last_temperature_update = fMCMCCurrentIteration;
    status.temperatures_frozen = status.temperature_updates == fPTTuning.rounds;
    WritePTTemperatureUpdate(old, rejection, update.damping);
    WritePTTemperatureConfiguration();
    BCLog::OutSummary(Form(" --> PT ladder update %u/%u at sweep %i%s; resetting proposal-training statistics.",
        status.temperature_updates, fPTTuning.rounds, fMCMCCurrentIteration, status.temperatures_frozen ? " (temperatures now frozen)" : ""));
    if (!status.temperatures_frozen) {
        fPTTuning.batch_sweeps *= 2;
        fPTTuning.next_sweep += fPTTuning.batch_sweeps;
    }
    std::fill(fPTTuning.rejection_sum.begin(), fPTTuning.rejection_sum.end(), 0.);
    std::fill(fPTTuning.observations.begin(), fPTTuning.observations.end(), 0);
    // Local proposals need to learn the new distributions. Retain their current
    // shapes, but restart the covariance gain and discard old-grid moments.
    fMultivariateCovarianceUpdates = 0;
    fMCMCNIterationsConvergenceGlobal = -1;
    for (auto& statistics : fMCMCStatistics) statistics.Reset(false, true);
    return true;
}

void BCEngineMCMC::WritePTTemperatureUpdate(const std::vector<double>& old,
                                           const std::vector<double>& rejections,
                                           double damping, bool initial) const
{
    if (fPTTuning.output.empty()) return;
    std::ofstream out(fPTTuning.output.c_str(), initial ? std::ios::trunc : std::ios::app);
    if (!out) throw std::runtime_error("Cannot open PT temperature history: " + fPTTuning.output);
    out.imbue(std::locale::classic()); out << std::setprecision(17);
    if (initial) out << "update,sweep,temperature_index,old_temperature,temperature,mean_rejection_to_next,observations,barrier,damping,frozen\n";
    double barrier = 0.; for (double r : rejections) barrier += r;
    for (unsigned t = 0; t < old.size(); ++t) {
        out << fWarmup.status.temperature_updates << ',' << fMCMCCurrentIteration << ',' << t << ','
            << old[t] << ',' << fPT.temperatures[t] << ',';
        if (t < rejections.size()) out << rejections[t];
        out << ',' << (t < fPTTuning.observations.size() ? fPTTuning.observations[t] : 0) << ','
            << barrier << ',' << damping << ',' << fWarmup.status.temperatures_frozen << '\n';
    }
    out.flush(); if (!out) throw std::runtime_error("Cannot flush PT temperature history: " + fPTTuning.output);
    out.close(); if (!out) throw std::runtime_error("Cannot close PT temperature history: " + fPTTuning.output);
}

void BCEngineMCMC::WritePTTemperatureConfiguration() const
{
    if (fPTDiagnostics.prefix.empty()) return;
    const auto path = fPTDiagnostics.prefix + "pt_temperatures.csv";
    std::ofstream out((path+".tmp").c_str());
    if (!out) throw std::runtime_error("Cannot write PT temperatures: " + path);
    out.imbue(std::locale::classic()); out << std::setprecision(17) << "slot,ladder,temperature\n";
    for (unsigned c = 0; c < GetNChains(); ++c)
        out << c << ',' << c%GetNColdChains() << ',' << GetChainTemperature(c) << '\n';
    out.flush(); if (!out) throw std::runtime_error("Cannot flush PT temperatures: " + path);
    out.close();
    if (!out || std::rename((path+".tmp").c_str(),path.c_str()) != 0)
        throw std::runtime_error("Cannot replace PT temperatures: " + path);
}
