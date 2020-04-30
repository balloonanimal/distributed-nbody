#include "simulation.h"
#include "bh_tree.h"
#include "mpi_utils.h"

Simulation *init_simulation(void) {
  Simulation *sim = calloc(1, sizeof(Simulation));
  if (sim == NULL) {
    perror("calloc");
    exit(1);
  }

  // Defaults
  sim->elapsed_time = 0.0;
  sim->dt = 0.1;
  sim->G = 1.0;
  sim->softening = 1E-6;
  sim->width = 10000;
  sim->gravity_method = PARTICLE_PARTICLE;
  sim->integration_method = LOCKSTEP;

  return sim;
}

void grow_array(ParticleArray *arr, int len) {
  while (len >= arr->allocation) {
    if (arr->allocation == 0) {
      arr->allocation = 128;
    } else {
      arr->allocation *= 2;
    }
    arr->particles =
        realloc(arr->particles, sizeof(Particle) * arr->allocation);
  }
}

void add_particle_to_array(ParticleArray *arr, Particle pt) {
  // grow array
  if (arr->len >= arr->allocation) {
    grow_array(arr, arr->len + 1);
  }
  arr->particles[arr->len] = pt;
  arr->len += 1;
}

void add_particle(Simulation *sim, Particle pt, bool owned) {
  add_particle_to_array(&sim->parray, pt);
  if (owned) {
    sim->N_owned += 1;
  }
}

void calc_gravity(Simulation *sim) {
  switch (sim->gravity_method) {
  case PARTICLE_PARTICLE:
    calc_gravity_pp(sim);
    break;
  case TREE:
    calc_gravity_tree(sim);
    break;
  default:
    printf("Unsupported gravity method\n");
    exit(1);
  }
}

void calc_gravity_pp(Simulation *sim) {
  Particle *particles = sim->parray.particles;
  int N = sim->parray.len;
  int N_owned = sim->N_owned;

  // zero all accelerations
  for (int i = 0; i < N_owned; i++) {
    particles[i].acceleration.x = 0;
    particles[i].acceleration.y = 0;
    particles[i].acceleration.z = 0;
  }

  // iterate through all pairs of particles
  // TODO: use plumber sphere method from AY190 paper
  //       the calculation for r is unneccesarily costly here

  // particles this process is responsible for updating
  for (int i = 0; i < N_owned; i++) {
    //  all particles in simulation
    for (int j = i + 1; j < N; j++) {
      double dx = particles[i].position.x - particles[j].position.x;
      double dy = particles[i].position.y - particles[j].position.y;
      double dz = particles[i].position.z - particles[j].position.z;
      double r2 = dx * dx + dy * dy + dz * dz + sim->softening * sim->softening;
      double almost_a = sim->G / r2;
      // update particle i
      double almost_ai = -almost_a * particles[j].mass;
      particles[i].acceleration.x = almost_ai * dx;
      particles[i].acceleration.y = almost_ai * dy;
      particles[i].acceleration.z = almost_ai * dz;

      if (j < N_owned) {
        // also update particle j
        double almost_aj = almost_a * particles[i].mass;
        particles[j].acceleration.x = almost_aj * dx;
        particles[j].acceleration.y = almost_aj * dy;
        particles[j].acceleration.z = almost_aj * dz;
      }
    }
  }
}

void calc_gravity_tree(Simulation *sim) {
  Particle *particles = sim->parray.particles;
  int N_owned = sim->N_owned;

  // zero all accelerations
  for (int i = 0; i < N_owned; i++) {
    particles[i].acceleration.x = 0;
    particles[i].acceleration.y = 0;
    particles[i].acceleration.z = 0;
  }

  // TODO: definitely doesn't belong here
  // HACK HACK HACK HACK HACK
  // FIXME FIXME FIXME FIXME
  free_tree(sim);
  build_tree(sim);

  // tree calculation
  Vec3d a;
  for (int i = 0; i < N_owned; i++) {
    // TODO maybe this sort of assignment is slower that doing it inside of the
    // fn?
    a = calc_gravity_tree_helper(sim, sim->tree, &particles[i]);
    particles[i].acceleration.x = a.x;
    particles[i].acceleration.y = a.y;
    particles[i].acceleration.z = a.z;
  }
}

Vec3d calc_gravity_tree_helper(Simulation *sim, BHTreeNode *node,
                               Particle *pt) {
  // node is terminal and contains pt, base case for recursion
  Vec3d a = {0, 0, 0};
  if (pt == node->pt) {
    return a;
  }

  double dx = (node->COM_position.x - pt->position.x);
  double dy = (node->COM_position.y - pt->position.y);
  double dz = (node->COM_position.z - pt->position.z);
  double r2 = dx * dx + dy * dy + dz * dz + sim->softening * sim->softening;
  double w2 = node->width * node->width;
  // Far away or lone node, estimate is good
  if (w2 / r2 < 1 || node->pcount == 1) {
    double almost_a = -sim->G * node->COM_mass / r2;
    a.x = almost_a * dx;
    a.y = almost_a * dy;
    a.z = almost_a * dz;
  }
  // Too close to estimate, recur to child nodes
  else {
    Vec3d a_;
    for (int oct = 0; oct < 8; oct++) {
      if (node->children[oct] != NULL) {
        a_ = calc_gravity_tree_helper(sim, node->children[oct], pt);
        a.x += a_.x;
        a.y += a_.y;
        a.z += a_.z;
      }
    }
  }
  return a;
}

void integrate(Simulation *sim, int time_steps) {
  for (int step = 1; step <= time_steps; step++) {
    integrate_step(sim);
  }
}

void integrate_step(Simulation *sim) {
  if (sim->use_mpi) {
    MPI_sync_particles(sim);
  }
  switch (sim->integration_method) {
  case LOCKSTEP:
    integrate_step_lock(sim);
    break;
  case LEAPFROG:
    integrate_step_leapfrog(sim);
    break;
  default:
    printf("Unsupported gravity method\n");
    exit(1);
  }
}

void integrate_step_lock(Simulation *sim) {
  Particle *particles = sim->parray.particles;
  int N_owned = sim->N_owned;
  double dt = sim->dt;

  // update accelerations
  calc_gravity(sim);

  // update velocities
  for (int i = 0; i < N_owned; i++) {
    particles[i].velocity.x += particles[i].acceleration.x * dt;
    particles[i].velocity.y += particles[i].acceleration.y * dt;
    particles[i].velocity.z += particles[i].acceleration.z * dt;
  }

  // update positions
  for (int i = 0; i < N_owned; i++) {
    particles[i].position.x += particles[i].velocity.x * dt;
    particles[i].position.y += particles[i].velocity.y * dt;
    particles[i].position.z += particles[i].velocity.z * dt;
  }

  sim->elapsed_time += dt;
}

void integrate_step_leapfrog(Simulation *sim) {
  int i;
  Particle *particles = sim->parray.particles;
  int N_owned = sim->N_owned;

  double dt = sim->dt;

  // first half update positions
  for (i = 0; i < N_owned; i++) {
    particles[i].position.x += 0.5 * particles[i].velocity.x * dt;
    particles[i].position.y += 0.5 * particles[i].velocity.y * dt;
    particles[i].position.z += 0.5 * particles[i].velocity.z * dt;
  }

  // update accelerations
  calc_gravity(sim);

  // update velocities
  for (i = 0; i < N_owned; i++) {
    particles[i].velocity.x += particles[i].acceleration.x * dt;
    particles[i].velocity.y += particles[i].acceleration.y * dt;
    particles[i].velocity.z += particles[i].acceleration.z * dt;
  }

  // second half update positions
  for (i = 0; i < N_owned; i++) {
    particles[i].position.x += 0.5 * particles[i].velocity.x * dt;
    particles[i].position.y += 0.5 * particles[i].velocity.y * dt;
    particles[i].position.z += 0.5 * particles[i].velocity.z * dt;
  }

  // TODO: MPI communicate updated positions / velocities

  sim->elapsed_time += dt;
}
