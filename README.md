# Distributed Simulation of the Nbody Problem

This repository contains the source code for my final project in Swarthmore's parallel and distributed computing class. In this project I worked to implement a program to simulate the nbody problem using two advanced approaches. The naive method requires the forces between all possible pairs of particles to be calculated, a total of n choose 2 or O(n^2) runtime. This fails to scale to the extrememly large systems of stars and other gravitational bodies that are relevant to modern astrophysics and cosmology, prompting both algorithmic and parallelization speedups to be developed.

## Barnes-Hut Tree
Instead of calculating the force on every pair of particles, the Barnes-Hut tree method recursively subdivides the simulation space in eighths, forming a tree based representation of the problem. Each node of the tree stores the total mass and center of mass for all the particles it contains. For particles sufficiently far away these summary values can be used to estimate the gravitational effect of all particles within the node. This method of grouping together the influence of far away particles results in a signifigant complexity reduction to O(n logn).

## Distributed Calculation
While the naive method listed earlier fails to scale well to large problem sizes, we can ameliorate this problem by splitting the work across several machines. For m machines, each machine will only have to calculate the gravitational forces on the n/m particles it is responsible for, reducing our O(n^2) loop to O(n^2/m). The dominating force in our current implementation is the communication overhead of sending updated particle positions at each timestep.

It is also possible to similarly parallelize the Barnes-Hut tree approach, which is the approach taken in actual research ready implementations of this problem, but this ended up being outside the scope of this project.

## Background
Some additional information can be found in the project proposal directory.
