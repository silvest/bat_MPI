#include "BCDisjointFlatPrior.h"
#include "BCAux.h"

#include <TH1.h>
#include <TRandom.h>

#include <algorithm>
#include <stdexcept>

BCDisjointFlatPrior::BCDisjointFlatPrior(std::vector<Interval> intervals)
    : fIntervals(std::move(intervals))
{
    if (fIntervals.empty())
        throw std::invalid_argument("A disjoint flat prior requires at least one interval");
    for (const auto& interval : fIntervals)
        if (!std::isfinite(interval.first) || !std::isfinite(interval.second) || interval.first >= interval.second)
            throw std::invalid_argument("Prior intervals require finite, increasing endpoints");

    std::sort(fIntervals.begin(), fIntervals.end());
    for (size_t i = 1; i < fIntervals.size(); ++i)
        if (fIntervals[i].first <= fIntervals[i-1].second)
            throw std::invalid_argument("Closed prior intervals must be disjoint");

    fTotalWidth = Width(fIntervals.front().first, fIntervals.back().second);
    fLogDensity = -std::log(fTotalWidth);
    fPriorFunction.SetRange(fIntervals.front().first, fIntervals.back().second);
}

BCDisjointFlatPrior::BCDisjointFlatPrior(const BCDisjointFlatPrior& other)
    : BCPrior(other),
      fIntervals(other.fIntervals),
      fTotalWidth(other.fTotalWidth),
      fLogDensity(other.fLogDensity)
{
    fPriorFunction.SetRange(other.fPriorFunction.GetXmin(), other.fPriorFunction.GetXmax());
}

double BCDisjointFlatPrior::GetLogPrior(double x)
{
    const auto it = std::lower_bound(fIntervals.begin(), fIntervals.end(), x,
        [](const Interval& interval, double value) { return interval.second < value; });
    return it != fIntervals.end() && it->first <= x ? fLogDensity : -INFINITY;
}

long double BCDisjointFlatPrior::Width(double xmin, double xmax) const
{
    long double width = 0;
    for (const auto& interval : fIntervals) {
        const long double lo = std::max(xmin, interval.first);
        const long double hi = std::min(xmax, interval.second);
        if (hi > lo)
            width += hi - lo;
    }
    return width;
}

double BCDisjointFlatPrior::GetIntegral(double xmin, double xmax)
{
    return Width(xmin, xmax) / fTotalWidth;
}

BCH1D BCDisjointFlatPrior::GetBCH1D(TH1* bins, const std::string& name)
{
    BCH1D result;
    if (!bins || bins->GetDimension() != 1)
        return result;

    std::unique_ptr<TH1> histogram(BCAux::OwnClone(bins, name));
    histogram->Reset();
    FillHistogramByIntegral(histogram.get());
    result.SetHistogram(histogram.get());
    result.SetLocalMode(GetMode(bins->GetXaxis()->GetXmin(), bins->GetXaxis()->GetXmax()));
    return result;
}

double BCDisjointFlatPrior::GetMode(double xmin, double xmax)
{
    for (const auto& interval : fIntervals) {
        const long double lo = std::max(xmin, interval.first);
        const long double hi = std::min(xmax, interval.second);
        if (lo <= hi)
            return (lo + hi) / 2;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

long double BCDisjointFlatPrior::Moment(unsigned n, double xmin, double xmax, long double center) const
{
    const long double width = Width(xmin, xmax);
    if (width == 0)
        return std::numeric_limits<long double>::quiet_NaN();
    if (n == 0)
        return 1;
    long double moment = 0;
    for (const auto& interval : fIntervals) {
        const long double lo = std::max(xmin, interval.first);
        const long double hi = std::min(xmax, interval.second);
        if (hi <= lo)
            continue;
        const long double a = lo - center;
        const long double b = hi - center;
        // The divided difference of x^(n+1) avoids subtracting nearby endpoint powers.
        long double sum = 1;
        long double power = 1;
        for (unsigned i = 0; i < n; ++i) {
            power *= a;
            sum = sum * b + power;
        }
        moment += ((hi - lo) / width) * sum / (n + 1.L);
    }
    return moment;
}

double BCDisjointFlatPrior::GetRawMoment(unsigned n, double xmin, double xmax)
{
    return Moment(n, xmin, xmax, 0);
}

double BCDisjointFlatPrior::GetCentralMoment(unsigned n, double xmin, double xmax)
{
    const long double mean = Moment(1, xmin, xmax, 0);
    if (n == 1 && std::isfinite(mean))
        return 0;
    return Moment(n, xmin, xmax, mean);
}

double BCDisjointFlatPrior::GetRandomValue(double xmin, double xmax, TRandom* const random)
{
    const long double width = Width(xmin, xmax);
    if (width == 0)
        return std::numeric_limits<double>::quiet_NaN();

    // Draw uniformly in allowed length, then insert the excluded gaps.
    long double offset = (random ? random : gRandom)->Rndm() * width;
    double last = std::numeric_limits<double>::quiet_NaN();
    for (const auto& interval : fIntervals) {
        const long double lo = std::max(xmin, interval.first);
        const long double hi = std::min(xmax, interval.second);
        if (hi <= lo)
            continue;
        if (offset < hi - lo)
            return lo + offset;
        offset -= hi - lo;
        last = hi;
    }
    return last; // Rounding at the upper endpoint.
}
