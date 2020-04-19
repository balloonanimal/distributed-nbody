#pragma once

#include "types.h"
#include "mpi_utils.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

// initialize simulation with default values
Simulation* init_simulation(void); 

// grow particle array to size len
void grow_array(ParticleArray* arr, int len);
// add particle to particlearray
void add_particle_to_array(ParticleArray* arr, Particle pt);
// add particle to simulation
void add_particle(Simulation* sim, Particle pt, bool owned);

// updates all particle accelerations according to the value of gravity_method
// TODO: does this statefully, maybe worth looking at immutability?
void calc_gravity(Simulation* sim);

// helper methods for individual gravity methods
void calc_gravity_pp(Simulation* sim);
void calc_gravity_tree(Simulation* sim);

// integrates the system for time_steps
// uses integrator based on value of integration_method
// calculates gravity based on balue of gravity_method
void integrate(Simulation* sim, int time_steps);

// calculate a single step of integration
void integrate_step(Simulation* sim);

// helper methods for integrators
void integrate_step_lock(Simulation* sim);
void integrate_step_leapfrog(Simulation* sim);

