// Communication-barrier inversion (Syed et al., Section 5), with bounded updates.
// Piecewise linear interpolation is monotone and needs no external spline library.
#ifndef BAT_PT_TEMPERATURE_SCHEDULE_H
#define BAT_PT_TEMPERATURE_SCHEDULE_H
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace BATPT {
struct ScheduleUpdate { std::vector<double> temperatures; double damping; };
inline ScheduleUpdate UpdateSchedule(const std::vector<double>& temperatures,
                                     const std::vector<double>& rejection, double damping)
{
    const unsigned n = temperatures.size();
    if (n < 3 || rejection.size()+1 != n || temperatures.front() != 1. ||
        !std::isfinite(damping) || damping <= 0. || damping > 1.)
        throw std::invalid_argument("Invalid PT temperature schedule update");
    std::vector<double> beta(n), mass(n, 0.), candidate(n);
    for (unsigned i = 0; i < n; ++i) {
        if (!std::isfinite(temperatures[i]) || temperatures[i] < 1. ||
            (i && temperatures[i] <= temperatures[i-1]))
            throw std::invalid_argument("Invalid PT temperature ladder");
        beta[n-1-i] = 1./temperatures[i];
    }
    for (unsigned i = 1; i < n; ++i) {
        const double r = rejection[n-1-i];
        if (!std::isfinite(r) || r < 0. || r > 1. || beta[i] <= beta[i-1])
            throw std::invalid_argument("Invalid PT rejection estimate or unresolved inverse-temperature gap");
        // Avoid singular inversion when all observed swaps would be accepted.
        mass[i] = mass[i-1] + std::max(1.e-12, r);
    }
    candidate.front() = beta.front(); candidate.back() = beta.back();
    for (unsigned i = 1, edge = 1; i+1 < n; ++i) {
        const double target = mass.back()*i/(n-1);
        while (edge+1 < n && mass[edge] < target) ++edge;
        candidate[i] = beta[edge-1] + (beta[edge]-beta[edge-1]) *
                       ((target-mass[edge-1])/(mass[edge]-mass[edge-1]));
    }
    // One common step length preserves endpoints and monotonicity, and bounds
    // EVERY gap, including the two endpoint gaps, rather than clipping knots.
    for (unsigned i = 1; i < n; ++i) {
        const double gap = beta[i]-beta[i-1], proposed = candidate[i]-candidate[i-1];
        if (proposed > gap) damping = std::min(damping, (std::exp(.5)-1.)*gap/(proposed-gap));
        if (proposed < gap) damping = std::min(damping, (1.-std::exp(-.5))*gap/(gap-proposed));
    }
    ScheduleUpdate result{temperatures, damping};
    for (unsigned i = 1; i+1 < n; ++i)
        result.temperatures[n-1-i] = 1./(beta[i] + damping*(candidate[i]-beta[i]));
    for (unsigned i = 1; i < n; ++i)
        if (!std::isfinite(result.temperatures[i]) || result.temperatures[i] <= result.temperatures[i-1])
            throw std::runtime_error("PT temperature update lost strict ordering");
    return result;
}
}
#endif
