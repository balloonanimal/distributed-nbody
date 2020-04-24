#include "simulation.h"
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

  // reproducible randomness
  srand(0);
  for (int i = 0; i < size; i++) {
    Particle p = {{0}};
    p.mass = rand() % 50;
    p.position.x = rand() % 100;
    p.position.y = rand() % 100;
    p.position.z = rand() % 100;
    add_particle(sim, p, true);
  }

  integrate(sim, steps);
  return 0;
}
