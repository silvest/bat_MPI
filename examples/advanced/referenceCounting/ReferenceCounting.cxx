#include "ReferenceCounting.h"

#include "RefPrior.h"
#include "GammaDistPrior.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>
#include <BAT/BCTF1Prior.h>
#include <BAT/BCTH1Prior.h>

#include <TH1.h>
#include <TF1.h>

// ---------------------------------------------------------
ReferenceCounting::ReferenceCounting(const std::string& name, unsigned nobs)
    : BCModel(name),
      fNObs(nobs)
{
    // define parameters
    AddParameter("s", 0, 40.); // index 0
    AddParameter("b", 0, 40.); // index 1
}

// ---------------------------------------------------------
bool ReferenceCounting::SetPrior(ReferenceCounting::EvalOption option, double shape, double rate)
{
  std::shared_ptr<RefPrior> refprior = std::make_shared<RefPrior>(shape, rate);
  std::shared_ptr<GammaDistPrior> conjprior = std::make_shared<GammaDistPrior>(shape, rate); // Poisson conjugate prior := Gamma Distribution

    // if analytic, use prior objects
    if (option == kAnalytic) {
        GetParameter("s").SetPrior(refprior);
        GetParameter("b").SetPrior(conjprior);
    } else {

        // both kHistogram and kApprox require signal histogram

        // create and fill signal prior
        TH1* h_s = GetParameter("s").CreateH1("hist_prior_s");
        refprior->FillHistogramByCenterValue(h_s);

        // don't need analytic refprior anymore
        refprior.reset();

        if (option == kHistogram) {

            // create and fill background prior
            TH1* h_b = GetParameter("b").CreateH1("hist_prior_b");
            conjprior->FillHistogramByCenterValue(h_b);
            // don't need analytic conjugate prior anymore
            conjprior.reset();

            GetParameter("s").SetPrior(std::make_shared<BCTH1Prior>(h_s));
            GetParameter("b").SetPrior(std::make_shared<BCTH1Prior>(h_b));

        } else if (option == kApprox) {

            // Create TF1 prior for signal and fit it to histogram created above
	  GetParameter("s").SetPrior(std::make_shared<BCTF1Prior>("sqrt(([0]*exp([1]*x^0.125))/(x+[0]*exp([1]*x^0.125)))",
                                       h_s->GetXaxis()->GetXmin(), h_s->GetXaxis()->GetXmax()));
            h_s->Fit(&GetParameter("s").GetPrior()->GetFunction());
            delete h_s;

            GetParameter("b").SetPrior(conjprior);

        } else {
	  conjprior.reset();
            return false;
        }
    }

    return GetParameter("s").GetPrior()->IsValid() && GetParameter("b").GetPrior()->IsValid();
}

// ---------------------------------------------------------
double ReferenceCounting::LogLikelihood(const std::vector<double>& parameters)
{
    // Likelihood = Poisson( nobs | expectation = s + b)
    return BCMath::LogPoisson(fNObs, parameters[0] + parameters[1]);
}
