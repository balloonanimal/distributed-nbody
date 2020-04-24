#include <stdio.h>

#include "mpi_utils.h"
#include "simulation.h"

int main(int argc, char *argv[]) {
  Simulation *sim = init_simulation();
  sim->gravity_method = PARTICLE_PARTICLE;
  sim->integration_method = LOCKSTEP;

  Particle pt1 = {{0}};
  pt1.mass = 100;
  add_particle(sim, pt1, true);

  Particle pt2 = {{0}};
  pt2.mass = 1;
  pt2.position.x = 100;
  add_particle(sim, pt2, true);

  sim->N_owned = 2;

  Particle *pt1_ = &(sim->parray.particles[0]);
  Particle *pt2_ = &(sim->parray.particles[1]);

  for (int i = 0; i < 5; i++) {
    integrate(sim, 1);
    printf("STEP %d\n", i);
    printf("pt1 position = %f\n", pt1_->position.x);
    printf("pt2 position = %f\n", pt2_->position.x);
  }

  /* printf("pt1 acceleration = x: %lf y: %lf z: %lf\n", */
  /*        pt1_->acceleration.x, pt1_->acceleration.y, pt1_->acceleration.z);
   */
  /* printf("pt2 acceleration = x: %lf y: %lf z: %lf\n", */
  /*        pt2_->acceleration.x, pt2_->acceleration.y, pt2_->acceleration.z);
   */
  /* // pt2 should fall in to pt1 */
  /* /\* integrate(sim, 1); *\/ */
  /* printf("Calculating gravity\n"); */
  /* calc_gravity(sim); */

  /* printf("pt1 acceleration = x: %lf y: %lf z: %lf\n", */
  /*        pt1_->acceleration.x, pt1_->acceleration.y, pt1_->acceleration.z);
   */
  /* printf("pt2 acceleration = x: %lf y: %lf z: %lf\n", */
  /*        pt2_->acceleration.x, pt2_->acceleration.y, pt2_->acceleration.z);
   */
  integrate(sim, 100);

  return 0;
}
