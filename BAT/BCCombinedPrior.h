#ifndef __BCCOMBINEDPRIOR__H
#define __BCCOMBINEDPRIOR__H

/*!
 * \class BCCombinedPrior
 * \brief A class to represent a combined prior with both Gaussian and Uniform uncertainties for a parameter
 * \author Luca Silvestrini
 * \version 1.0
 * \date 05.2024
 * \ingroup Priors
 */


/*
 * Copyright (C) 2024, the HEPfit core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCPrior.h"

#include <limits>
#include <cmath>

class TF1;

// ---------------------------------------------------------

class BCCombinedPrior : public BCPrior
{
public:
    /** \name Constructor & Destructor */
    /** @{ **/

    /** Constructor */
    BCCombinedPrior(double mean, double errGauss, double errUniform);

    /** Destructor */
    virtual ~BCCombinedPrior() {};

    /** @} **/

    /** \name Functions overloaded from BCPrior **/
    /** @{ **/

    /** Clone function */
    virtual BCPrior* Clone() const
    { return new BCCombinedPrior(*this); }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const
    { return std::isfinite(fMean) and std::isfinite(fSigma) and fSigma > 0 and std::isfinite(fErrUniform) and fErrUniform > 0 and fPriorFunction.IsValid();}

    /**
     * Get log of prior
     * @param x value to evaluate log of prior at
     * @return log of prior */
    virtual double GetLogPrior(double x)
    { return log(fPriorFunction.Eval(x)); }

    /**
     * Get prior
     * @param x value to evaluate prior at
     * @param normalize Whether to normalize prior with stored integral
     * @return prior */
    virtual double GetPrior(double x, bool normalize = false)
    { return fPriorFunction.Eval(x) * ((normalize) ? exp(-fLogIntegral) : 1); }

    /**
     * Return mode of prior (in range).
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mode of prior in range. */
    virtual double GetMode(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { if (fMean < xmin) return xmin; if (fMean > xmax) return xmax; return fMean; }

    /** @} **/

    /** \name Setters **/
    /** @{ **/

    void SetMean(double mean)
    { fMean = mean; }

    void SetErrGauss(double errGauss)
    { fSigma = errGauss; }

    void SetErrUniform(double errUniform)
    { fErrUniform = errUniform; }

    void SetParameters(double mean, double errGauss, double errUniform)
    { SetMean(mean); SetErrGauss(errGauss); SetErrUniform(errUniform); fMax = fMean + 5*fSigma + fErrUniform; fMin = fMean - 5*fSigma - fErrUniform;}

    void SetFunction();

    /** @} **/

protected:
    double fMean=0.;									///< mean of Gaussian
    double fSigma=0.;								///< std dev of Gaussian
    double fErrUniform=0.;						///< error of uniform prior
    double fMin=0.;									///< minimum value of prior
    double fMax=0.;									///< maximum value of prior

};

#endif
