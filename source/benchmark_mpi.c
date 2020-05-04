#include "simulation.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    printf("Usage: ./benchmark [num_particles] [steps]\n");
    exit(1);
  }
  int size = atoi(argv[1]);
  int steps = atoi(argv[2]);

  Simulation *sim = init_simulation();
  sim->gravity_method = PARTICLE_PARTICLE;
  sim->integration_method = LOCKSTEP;

  // start MPI
  MPI_setup(sim);

  if (size % sim->MPI_pcount != 0) {
    printf("we only currently support sizes that are divided evenly by the "
           "number of processes\n");
    exit(1);
  }

  int split_size = size / sim->MPI_pcount;
  // reproducible randomness
  srand(100 * sim->MPI_rank);
  for (int i = 0; i < split_size; i++) {
    Particle p = {0};
    p.mass = rand() % 50;
    p.pos_x = rand() % (int)(2 * sim->width) - sim->width;
    p.pos_y = rand() % (int)(2 * sim->width) - sim->width;
    p.pos_z = rand() % (int)(2 * sim->width) - sim->width;
    add_particle(sim, p, true);
  }
  print_array(sim);

  integrate(sim, steps);

  print_array(sim);

  MPI_teardown(sim);
  return 0;
}
