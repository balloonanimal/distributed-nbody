#include <mpi.h>
#include <stddef.h>
#include "types.h"

void MPI_setup(Simulation* sim);
void MPI_teardown(Simulation* sim);
void MPI_sync_particles(Simulation* sim);
