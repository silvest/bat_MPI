// Bounded proposal adaptation followed by discarded fixed-kernel settling.
#include "BCEngineMCMC.h"
#include <cstdio>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <limits>
#include <locale>
#include <stdexcept>

void BCEngineMCMC::SetPreRunPolicy(PreRunPolicy policy, unsigned settling_sweeps)
{
    if (fMCMCPhase != kUnsetPhase || !fMCMCStates.empty())
        throw std::logic_error("Configure the prerun policy before MCMC initialization");
    if ((policy != kPreRunLegacy && policy != kPreRunStaged) ||
        settling_sweeps > static_cast<unsigned>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Invalid prerun policy or settling sweep count");
    fWarmup.policy = policy;
    fWarmup.settling_sweeps = settling_sweeps;
    ResetPreRunStatus();
}

void BCEngineMCMC::SetScaleAdaptation(ScaleAdaptation method, double target, double gain, double decay)
{
    if (fMCMCPhase != kUnsetPhase || !fMCMCStates.empty())
        throw std::logic_error("Configure scale adaptation before MCMC initialization");
    if ((method != kScaleLegacy && method != kScaleSmooth) || !std::isfinite(target) || target <= 0. || target >= 1. ||
        !std::isfinite(gain) || gain <= 0. || !std::isfinite(decay) || decay <= 0.5 || decay > 1.)
        throw std::invalid_argument("Invalid scale adaptation method, target, gain or decay");
    fWarmup.scale_adaptation = method;
    fWarmup.scale_target = target; fWarmup.scale_gain = gain; fWarmup.scale_decay = decay;
    ResetPreRunStatus();
}

BCEngineMCMC::ScaleUpdateResult BCEngineMCMC::UpdateSmoothProposalScales(bool joint)
{
    if (fMCMCPhase != kPreRun || fWarmup.fixed_kernel || fWarmup.status.stage != "adapting")
        throw std::logic_error("Scale adaptation is allowed only during adaptive warmup");
    auto& updates = joint ? fWarmup.status.joint_scale_updates : fWarmup.status.factorized_scale_updates;
    const double eta = fWarmup.scale_gain / std::pow(double(++updates), fWarmup.scale_decay);
    const double lower = std::max(fMCMCScaleFactorLowerLimit, std::numeric_limits<double>::min());
    const double upper = fMCMCScaleFactorUpperLimit;
    ScaleUpdateResult result;
    for (unsigned c = 0; c < GetNChains(); ++c) {
        if (!fMCMCStatistics[c].n_samples_efficiency) continue;
        for (unsigned p = 0; p < (joint ? 1 : GetNParameters()); ++p) {
            if (!joint && GetParameter(p).Fixed()) continue;
            const double rate = fMCMCStatistics[c].efficiency[p];
            double& scale = fMCMCProposalFunctionScaleFactor[c][p];
            if (!std::isfinite(rate) || rate < 0. || rate > 1. || !std::isfinite(scale) || scale <= 0.)
                throw std::runtime_error("Smooth adaptation requires finite positive scales and finite acceptance estimates");
            const double old = scale;
            const bool outlier = rate < fMCMCEfficiencyMin || rate > fMCMCEfficiencyMax;
            result.all_efficient &= !outlier;
            if (outlier && ((rate < fWarmup.scale_target && scale > lower) ||
                            (rate > fWarmup.scale_target && scale < upper))) result.adjustable = true;
            // Work in log standard deviation. The joint native multiplier is a
            // variance factor, unlike the coordinate proposal's standard deviation.
            const double delta = std::max(-0.5, std::min(0.5, eta * (rate - fWarmup.scale_target)));
            const double log_scale = std::max(std::log(lower),
                std::min(std::log(upper), std::log(scale) + (joint ? 2. : 1.) * delta));
            scale = std::max(lower, std::min(upper, std::exp(log_scale)));
            result.changed |= scale != old;
        }
    }
    return result;
}

void BCEngineMCMC::SetPreRunStatusOutput(const std::string& filename)
{
    if (fMCMCPhase != kUnsetPhase || !fMCMCStates.empty())
        throw std::logic_error("Configure warmup status output before MCMC initialization");
    fWarmup.output = filename;
}

void BCEngineMCMC::ResetPreRunStatus()
{
    fWarmup.status = PreRunStatus();
    fWarmup.status.policy = fWarmup.policy;
    fWarmup.status.scale_adaptation = fWarmup.scale_adaptation;
    fWarmup.status.scale_target = fWarmup.scale_target;
    fWarmup.status.scale_gain = fWarmup.scale_gain;
    fWarmup.status.scale_decay = fWarmup.scale_decay;
    fWarmup.fixed_kernel = false;
}

void BCEngineMCMC::UpdatePreRunStatusOutput() const
{
    if (!fWarmup.output.empty()) WritePreRunStatus(fWarmup.output);
}

void BCEngineMCMC::WritePreRunStatus(const std::string& filename) const
{
    if (filename.empty()) throw std::invalid_argument("Warmup status requires a filename");
    const std::string temporary = filename + ".tmp";
    std::ofstream out(temporary.c_str(), std::ios::trunc);
    if (!out) throw std::runtime_error("Cannot open warmup status: " + temporary);
    out.imbue(std::locale::classic());
    out << std::setprecision(17);
    const auto& s = fWarmup.status;
    out << "policy,stage,stop_reason,adaptation_sweeps,settling_sweeps,adaptation_complete,settling_complete,proposal_loaded,"
           "adaptation_efficiency_outliers,adaptation_zero_acceptance,settling_efficiency_outliers,settling_zero_acceptance,"
           "legacy_rvalue_checked,legacy_rvalue_passed,posterior_convergence_assessed,scale_adaptation,scale_target,scale_gain,scale_decay,scale_max_log_step,factorized_scale_updates,joint_scale_updates\n"
        << (s.policy == kPreRunStaged ? "staged" : "legacy") << ',' << s.stage << ',' << s.stop_reason << ','
        << s.adaptation_sweeps << ',' << s.settling_sweeps << ',' << s.adaptation_complete << ','
        << s.settling_complete << ',' << s.proposal_loaded << ',' << s.adaptation_efficiency_outliers << ','
        << s.adaptation_zero_acceptance << ',' << s.settling_efficiency_outliers << ','
        << s.settling_zero_acceptance << ',' << s.legacy_rvalue_checked << ',' << s.legacy_rvalue_passed << ','
        << s.posterior_convergence_assessed << ',' << (s.scale_adaptation == kScaleSmooth ? "smooth" : "legacy") << ','
        << s.scale_target << ',' << s.scale_gain << ',' << s.scale_decay << ",0.5,"
        << s.factorized_scale_updates << ',' << s.joint_scale_updates << '\n';
    out.flush();
    if (!out) throw std::runtime_error("Cannot flush warmup status: " + temporary);
    out.close();
    if (!out || std::rename(temporary.c_str(), filename.c_str()) != 0)
        throw std::runtime_error("Cannot replace warmup status: " + filename);
}

void BCEngineMCMC::RecordPreRunEfficiencies(bool settling, bool joint)
{
    unsigned outliers = 0, zero = 0;
    for (unsigned c = 0; c < GetNChains(); ++c) {
        if (!fMCMCStatistics[c].n_samples_efficiency) continue;
        for (unsigned p = 0; p < (joint ? 1 : GetNParameters()); ++p) {
            if (!joint && GetParameter(p).Fixed()) continue;
            const double rate = fMCMCStatistics[c].efficiency[p];
            if (!std::isfinite(rate) || rate < fMCMCEfficiencyMin || rate > fMCMCEfficiencyMax) ++outliers;
            if (rate == 0.) ++zero;
        }
    }
    if (settling) {
        fWarmup.status.settling_efficiency_outliers = outliers;
        fWarmup.status.settling_zero_acceptance = zero;
    } else {
        fWarmup.status.adaptation_efficiency_outliers = outliers;
        fWarmup.status.adaptation_zero_acceptance = zero;
    }
}

void BCEngineMCMC::RunFixedProposalWarmup(unsigned sweeps)
{
    fWarmup.status.stage = "settling";
    fWarmup.status.settling_sweeps = 0;
    fWarmup.status.settling_complete = false;
    // A separate phase uses fresh efficiency counters and journeys. No adaptation,
    // user iteration callbacks, histograms, ROOT draws or convergence statistics.
    for (auto& statistics : fMCMCStatistics) statistics.Reset(false, true);
    ResetReplicaExchangeStatistics();
    fMCMCPhase = kPreRun;
    fMCMCCurrentIteration = 0;
    fWarmup.fixed_kernel = true;
    struct Restore { bool& flag; ~Restore() { flag = false; } } restore{fWarmup.fixed_kernel};
    UpdatePreRunStatusOutput();
    BeginPTDiagnostics("warmup");
    BCLog::OutSummary(Form(" --> Discarding %u settling sweeps with fixed proposals; posterior convergence must be assessed separately.", sweeps));
    while (fMCMCCurrentIteration < static_cast<int>(sweeps)) {
        GetNewPointMetropolis();
        fWarmup.status.settling_sweeps = fMCMCCurrentIteration;
        if (fMCMCNIterationsPreRunCheck && fMCMCCurrentIteration % fMCMCNIterationsPreRunCheck == 0) {
            RecordPreRunEfficiencies(true, UseMultivariateProposal());
            UpdatePreRunStatusOutput();
        }
    }
    RecordPreRunEfficiencies(true, UseMultivariateProposal());
    EndPTDiagnostics();
    fWarmup.status.settling_complete = true;
    fWarmup.status.stage = "complete";
    if (fWarmup.status.adaptation_efficiency_outliers || fWarmup.status.settling_efficiency_outliers)
        BCLog::OutWarning(Form("Warmup proposal diagnostics: %u acceptance estimates outside target in the last adaptation window; %u during fixed settling. These are tuning diagnostics, not posterior convergence tests.",
            fWarmup.status.adaptation_efficiency_outliers, fWarmup.status.settling_efficiency_outliers));
    if (fWarmup.status.settling_zero_acceptance)
        BCLog::OutWarning(Form("Warmup proposal problems: %u parameter/slot kernels accepted no local moves during settling. Inspect the proposal and transport diagnostics before using production results.",
            fWarmup.status.settling_zero_acceptance));
    // Discard settling efficiencies as well as samples; production starts clean.
    for (auto& statistics : fMCMCStatistics) statistics.Reset(false, true);
    UpdatePreRunStatusOutput();
}
