#include "BCCombinedPrior.h"
#include "TF1.h"

BCCombinedPrior::BCCombinedPrior(double mean, double errGauss, double errUniform)
    : fMean(mean), fSigma(errGauss), fErrUniform(errUniform) 
{
    fMax = fMean + 5*fSigma + fErrUniform;
    fMin = fMean - 5*fSigma - fErrUniform;
    SetFunction();
}

void BCCombinedPrior::SetFunction()
    { 
        fPriorFunction = TF1("fPriorFunction", "(TMath::Erf((x-[0]+[2])/sqrt(2.)/[1])-TMath::Erf((x-[0]-[2])/sqrt(2.)/[1]))/4./[2]", fMin, fMax);
        fPriorFunction.SetParameter(0, fMean);
        fPriorFunction.SetParameter(1, fSigma);
        fPriorFunction.SetParameter(2, fErrUniform);
        fLogIntegral = log(fPriorFunction.Integral(fMin, fMax)); 
    }
