#ifndef BCDISJOINTFLATPRIOR_H
#define BCDISJOINTFLATPRIOR_H

#include "BCPrior.h"

#include <utility>
#include <vector>

/** \ingroup Priors
 * Uniform probability density on a union of disjoint, closed intervals.
 * Each interval's probability is its width divided by the total allowed width.
 */
class BCDisjointFlatPrior : public BCPrior
{
public:
    using Interval = std::pair<double, double>;

    /** Accept intervals in any order. Throws std::invalid_argument for an empty
     * union, nonfinite bounds, nonpositive widths, touching or overlapping intervals. */
    explicit BCDisjointFlatPrior(std::vector<Interval> intervals);
    BCDisjointFlatPrior(const BCDisjointFlatPrior& other);
    BCDisjointFlatPrior& operator=(const BCDisjointFlatPrior&) = delete;

    BCDisjointFlatPrior* Clone() const { return new BCDisjointFlatPrior(*this); }
    const std::vector<Interval>& GetIntervals() const { return fIntervals; }
    bool IsValid() const override { return true; }
    double GetLogPrior(double x) override;

    double GetIntegral(double xmin = -INFINITY, double xmax = INFINITY) override;
    double GetMode(double xmin = -INFINITY, double xmax = INFINITY) override;
    double GetRawMoment(unsigned n, double xmin = -INFINITY, double xmax = INFINITY) override;
    double GetCentralMoment(unsigned n, double xmin = -INFINITY, double xmax = INFINITY) override;
    double GetRandomValue(double xmin, double xmax, TRandom* const random = nullptr) override;
    BCH1D GetBCH1D(TH1* bins, const std::string& name = "prior") override;

private:
    long double Width(double xmin, double xmax) const;
    long double Moment(unsigned n, double xmin, double xmax, long double center) const;

    std::vector<Interval> fIntervals;
    long double fTotalWidth;
    double fLogDensity;
};

#endif
