#pragma once

#include <mpi.h>
#include <stdbool.h>
#include <stdint.h>

typedef struct {
  double x;
  double y;
  double z;
} Vec3d;

typedef struct {
  Vec3d position;
  Vec3d velocity;
  Vec3d acceleration;
  double mass;
} Particle;

typedef struct {
  Particle *particles;
  int len;
  int allocation;
} ParticleArray;

typedef enum { PARTICLE_PARTICLE = 0, TREE = 1 } GravityMethod;

typedef enum { LOCKSTEP = 0, LEAPFROG = 1 } IntegrationMethod;

typedef struct {
  double elapsed_time;  // total elapsed time simulation has run
  double dt;            // size of timestep
  double G;             // value of gravitational constant
  double softening;     // value of softening constant
  ParticleArray parray; // array of particles
  uint64_t N_owned;     // total number of particles this node owns
  GravityMethod
      gravity_method; // which algorithm to use for gravity calculations
  IntegrationMethod
      integration_method; // which algorithm to use for integration
  // TODO a bunch of MPI stuff
  bool use_mpi;
  int MPI_pcount;
  int MPI_rank;
  MPI_Datatype MPI_particle_type;
  ParticleArray *send_arrays;
  ParticleArray *recv_arrays;
  int *recv_sizes;
} Simulation;
