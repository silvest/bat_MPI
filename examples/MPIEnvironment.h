#ifndef BAT_EXAMPLE_MPI_ENVIRONMENT_H
#define BAT_EXAMPLE_MPI_ENVIRONMENT_H
#include <mpi.h>
// Examples have no distributed likelihood worker loop: launch one MPI process.
struct BATExampleMPIEnvironment {
    BATExampleMPIEnvironment() { MPI_Init(nullptr, nullptr); }
    ~BATExampleMPIEnvironment() { MPI_Finalize(); }
};
#endif
