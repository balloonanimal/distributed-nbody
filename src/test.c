#include <stdio.h>

#include "simulation.h"

int main(int argc, char *argv[]) {
    Simulation* sim = init_simulation();
    sim->gravity_method= PARTICLE_PARTICLE;
    sim->integration_method = LOCKSTEP;

    Particle pt1 = {};
    pt1.mass = 100;
    add_particle(sim, pt1);

    Particle pt2 = {};
    pt2.position.x = 100;
    /* pt2.velocity.y = 100; */
    pt2.mass = 1;
    add_particle(sim, pt1);

    // pt2 should fall in to pt1
    integrate(sim, 1);

    pt1 = sim->particles[0];
    pt2 = sim->particles[1];
    printf("pt1 position = x: %lf y: %lf z: %lf\n",
           pt1.position.x, pt1.position.y, pt1.position.z);
    printf("pt2 position = x: %lf y: %lf z: %lf\n",
           pt2.position.x, pt2.position.y, pt2.position.z);
    /* integrate(sim, 100); */


    return 0;
}
