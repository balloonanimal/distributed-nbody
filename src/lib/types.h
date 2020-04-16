#pragma once

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

typedef enum {
    PARTICLE_PARTICLE = 0,
    TREE = 1
} GravityMethod;

typedef enum {
    LOCKSTEP = 0,
    LEAPFROG = 1
} IntegrationMethod;

typedef struct {
    double elapsed_time; // total elapsed time simulation has run
    double dt; // size of timestep
    double G; // value of gravitational constant
    double softening; // value of softening constant
    Particle* particles; // array of particles
    uint64_t N; // total number of particles
    uint64_t particles_allocation;
    GravityMethod gravity_method; // which algorithm to use for gravity calculations
    IntegrationMethod integration_method; // which algorithm to use for integration
    // TODO a bunch of MPI stuff
} Simulation;
