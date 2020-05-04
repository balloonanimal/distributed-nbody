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
  sim->G = 1;
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
    particles[i].acc_x = 0;
    particles[i].acc_y = 0;
    particles[i].acc_z = 0;
  }

  // iterate through all pairs of particles
  // TODO: use plumber sphere method from AY190 paper
  //       the calculation for r is unneccesarily costly here

  // particles this process is responsible for updating
  for (int i = 0; i < N_owned; i++) {
    //  all particles in simulation
    for (int j = i + 1; j < N; j++) {
      double dx = particles[i].pos_x - particles[j].pos_x;
      double dy = particles[i].pos_y - particles[j].pos_y;
      double dz = particles[i].pos_z - particles[j].pos_z;
      double r2 = dx * dx + dy * dy + dz * dz + sim->softening * sim->softening;
      double almost_a = sim->G / r2;
      // update particle i
      double almost_ai = -almost_a * particles[j].mass;
      particles[i].acc_x = almost_ai * dx;
      particles[i].acc_y = almost_ai * dy;
      particles[i].acc_z = almost_ai * dz;

      if (j < N_owned) {
        // also update particle j
        double almost_aj = almost_a * particles[i].mass;
        particles[j].acc_x = almost_aj * dx;
        particles[j].acc_y = almost_aj * dy;
        particles[j].acc_z = almost_aj * dz;
      }
    }
  }
}

void calc_gravity_tree(Simulation *sim) {
  Particle *particles = sim->parray.particles;
  int N_owned = sim->N_owned;

  // zero all accelerations
  for (int i = 0; i < N_owned; i++) {
    particles[i].acc_x = 0;
    particles[i].acc_y = 0;
    particles[i].acc_z = 0;
  }

  // build tree
  free_tree(sim);
  build_tree(sim);

  // tree calculation
  Vec3d a;
  for (int i = 0; i < N_owned; i++) {
    a = calc_gravity_tree_helper(sim, sim->tree, &particles[i]);
    particles[i].acc_x = a.x;
    particles[i].acc_y = a.y;
    particles[i].acc_z = a.z;
  }
}

Vec3d calc_gravity_tree_helper(Simulation *sim, BHTreeNode *node,
                               Particle *pt) {
  // node is terminal and contains pt, base case for recursion
  Vec3d a = {0, 0, 0};
  if (pt == node->pt) {
    return a;
  }

  double dx = (node->COM_x - pt->pos_x);
  double dy = (node->COM_y - pt->pos_y);
  double dz = (node->COM_z - pt->pos_z);
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
    particles[i].vel_x += particles[i].acc_x * dt;
    particles[i].vel_y += particles[i].acc_y * dt;
    particles[i].vel_z += particles[i].acc_z * dt;
  }

  // update positions
  for (int i = 0; i < N_owned; i++) {
    particles[i].pos_x += particles[i].vel_x * dt;
    particles[i].pos_y += particles[i].vel_y * dt;
    particles[i].pos_z += particles[i].vel_z * dt;
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
    particles[i].pos_x += 0.5 * particles[i].vel_x * dt;
    particles[i].pos_y += 0.5 * particles[i].vel_y * dt;
    particles[i].pos_z += 0.5 * particles[i].vel_z * dt;
  }

  // update accelerations
  calc_gravity(sim);

  // update velocities
  for (i = 0; i < N_owned; i++) {
    particles[i].vel_x += particles[i].acc_x * dt;
    particles[i].vel_y += particles[i].acc_y * dt;
    particles[i].vel_z += particles[i].acc_z * dt;
  }

  // second half update positions
  for (i = 0; i < N_owned; i++) {
    particles[i].pos_x += 0.5 * particles[i].vel_x * dt;
    particles[i].pos_y += 0.5 * particles[i].vel_y * dt;
    particles[i].pos_z += 0.5 * particles[i].vel_z * dt;
  }

  sim->elapsed_time += dt;
}

// UTILS
void print_array(Simulation *sim) {
  int len = sim->parray.len;
  Particle *particles = sim->parray.particles;
  printf("%d: [", sim->MPI_rank);
  if (len == 0) {
    printf("]\n");
    return;
  }
  for (int i = 0; i < len; i++) {
    printf("(%.3f, %.3f, %.3f : %f)", particles[i].pos_x, particles[i].pos_y,
           particles[i].pos_z, particles[i].mass);
    if (i == sim->N_owned - 1) {
      printf(" | ");
    }
    if (i == len - 1) {
      printf("]\n");
    } else {
      printf(" ");
    }
  }
}
