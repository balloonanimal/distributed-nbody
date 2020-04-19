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

  // TODO: probably want to spin this out into its own fn
  ParticleArray gen_parray = {0};
  if (sim->MPI_rank == 0) {
    // reproducible randomness
    srand(0);
    for (int i = 0; i < size; i++) {
      Particle p = {0};
      p.mass = rand() % 50;
      p.position.x = rand() % 100;
      p.position.y = rand() % 100;
      p.position.z = rand() % 100;
      add_particle_to_array(&gen_parray, p);
    }
  }

  // grow array to fit particles
  int split_size = size / sim->MPI_pcount;
  grow_array(&sim->parray, split_size);

  // scatter particles to all processes
  MPI_Scatter(gen_parray.particles, split_size, sim->MPI_particle_type,
              sim->parray.particles, split_size, sim->MPI_particle_type, 0,
              MPI_COMM_WORLD);

  integrate(sim, steps);
  MPI_teardown(sim);
  return 0;
}
