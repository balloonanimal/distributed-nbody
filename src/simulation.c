#include "simulation.h"

Simulation* init_simulation(void) {
    Simulation *sim = calloc(0, sizeof(Simulation));

    // Defaults
    sim->elapsed_time = 0;
    sim->dt = 0.1;
    sim->G = 1;
    sim->softening = 1E-8;
    sim->N = 0;
    sim->particles_allocation = 0;
    sim->gravity_method = PARTICLE_PARTICLE;
    sim->integration_method = LOCKSTEP;

    return sim;
}

void add_particle(Simulation* sim, Particle pt) {
    // grow array
    if (sim->N >= sim->particles_allocation) {
        if (sim->particles_allocation == 0) {
            sim->particles_allocation = 128;
        } else {
            sim->particles_allocation *= 2;
        }
        sim->particles = realloc(sim->particles,
                                 sizeof(Particle) * sim->particles_allocation);
    }
    sim->particles[sim->N + 1] = pt;
    sim->N += 1;
}

void calc_gravity(Simulation* sim) {
    switch (sim->gravity_method) {
        case PARTICLE_PARTICLE:
            calc_gravity_pp(sim);
            break;
        case TREE:
            printf("Tree calculation not yet implemented\n");
            exit(1);
        default:
            printf("Unsupported gravity method\n");
            exit(1);
    }
}

void calc_gravity_pp(Simulation* sim) {
    int i, j;
    Particle* particles = sim->particles;

    // zero all accelerations
    for (i = 0; i < sim->N; i++){
        particles[i].acceleration.x = 0;
        particles[i].acceleration.y = 0;
        particles[i].acceleration.z = 0;
    }

    // iterate through all pairs of particles
    // TODO: use plumber sphere method from AY190 paper
    //       the calculation for r is unneccesarily costly here
    for (i = 0; i < sim->N - 1; i++) {
        for (j = i + 1; j < sim->N; j++) {
            double dx = particles[i].position.x - particles[j].position.x;
            double dy = particles[i].position.y - particles[j].position.y;
            double dz = particles[i].position.z - particles[j].position.z;
            double r = sqrt(dx*dx + dy*dy + dz*dz) + sim->softening;
            double almost_a = sim->G / (r*r);
            double almost_ai = almost_a * particles[j].mass;
            double almost_aj = -almost_a * particles[i].mass;

            particles[i].acceleration.x = almost_ai * dx;
            particles[i].acceleration.y = almost_ai * dy;
            particles[i].acceleration.z = almost_ai * dz;
            particles[j].acceleration.x = almost_aj * dx;
            particles[j].acceleration.y = almost_aj * dy;
            particles[j].acceleration.z = almost_aj * dz;
        }
    }
}

// TODO: implement
void calc_gravity_tree(Simulation* sim) {

}

void integrate(Simulation* sim, int time_steps) {
    int step;
    for(step = 1; step <= time_steps; step++) {
        integrate_step(sim);
    }
}

void integrate_step(Simulation* sim) {
    switch (sim->integration_method) {
        case LOCKSTEP:
            integrate_step_lock(sim);
            break;
        case LEAPFROG:
            printf("Leapfrog integration not yet implemented\n");
            exit(1);
        default:
            printf("Unsupported gravity method\n");
            exit(1);
    }
}

void integrate_step_lock(Simulation* sim) {
    int i;
    Particle* particles = sim->particles;
    double dt = sim->dt;

    // update accelerations
    calc_gravity(sim);

    // update velocities
    for (i = 0; i < sim->N; i++) {
        particles[i].velocity.x += particles[i].acceleration.x * dt;
        particles[i].velocity.y += particles[i].acceleration.y * dt;
        particles[i].velocity.z += particles[i].acceleration.z * dt;
    }

    // update positions
    for (i = 0; i < sim->N; i++) {
        particles[i].position.x += particles[i].velocity.x * dt;
        particles[i].position.y += particles[i].velocity.y * dt;
        particles[i].position.z += particles[i].velocity.z * dt;
    }

    sim->elapsed_time += dt;
}

void integrate_step_leapfrog(Simulation* sim) {
    int i;
    Particle* particles = sim->particles;
    double dt = sim->dt;

    // first half update positions
    for (i = 0; i < sim->N; i++) {
        particles[i].position.x += 0.5 * particles[i].velocity.x * dt;
        particles[i].position.y += 0.5 * particles[i].velocity.y * dt;
        particles[i].position.z += 0.5 * particles[i].velocity.z * dt;
    }

    // update accelerations
    calc_gravity(sim);

    // update velocities
    for (i = 0; i < sim->N; i++) {
        particles[i].velocity.x += particles[i].acceleration.x * dt;
        particles[i].velocity.y += particles[i].acceleration.y * dt;
        particles[i].velocity.z += particles[i].acceleration.z * dt;
    }

    // second half update positions
    for (i = 0; i < sim->N; i++) {
        particles[i].position.x += 0.5 * particles[i].velocity.x * dt;
        particles[i].position.y += 0.5 * particles[i].velocity.y * dt;
        particles[i].position.z += 0.5 * particles[i].velocity.z * dt;
    }

    sim->elapsed_time += dt;
}
