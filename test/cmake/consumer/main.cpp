#include <BAT/BCModel.h>
#include <BAT/BCParameter.h>
#include <mpi.h>
#ifdef TEST_MODELS
#include <BAT/BCGraphFitter.h>
#include <TGraphErrors.h>
#include <TF1.h>
#endif
#ifdef TEST_MTF
#include <BAT/BCMTF.h>
#endif
#include <stdexcept>
struct InstalledModel : BCModel {
    InstalledModel() : BCModel("installed_consumer") {}
    double LogLikelihood(const std::vector<double>& x) override { return -x.at(0)*x.at(0); }
};
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    {
        InstalledModel model;
        model.AddParameter("x", -2., 2.);
        model.SetPriorConstantAll();
        model.SetParallelTempering({1., 2.}, 2, 5);
        model.SetPTDiagnosticsOutput("consumer_stream_",10,0);
        if (model.GetNChains() != 4 || model.GetNColdChains() != 2)
            throw std::runtime_error("Installed PT API mismatch");
#ifdef TEST_MODELS
        TGraphErrors graph(2);
        graph.SetPoint(0, 0., 0.);
        graph.SetPoint(1, 1., 1.);
        TF1 function("line", "[0]*x", 0., 1.);
        BCGraphFitter fitter(graph, function);
#endif
#ifdef TEST_MTF
        BCMTF templates;
#endif
    }
    MPI_Finalize();
}
