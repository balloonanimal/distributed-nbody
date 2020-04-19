#include "simulation.h"
#include "types.h"

void MPI_setup(Simulation *sim) {
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &sim->MPI_pcount);
  MPI_Comm_rank(MPI_COMM_WORLD, &sim->MPI_rank);

  // Init buffers for all other processes
  sim->send_arrays = calloc(sim->MPI_pcount, sizeof(ParticleArray));
  sim->recv_arrays = calloc(sim->MPI_pcount, sizeof(ParticleArray));
  for (int pnum = 0; pnum < sim->MPI_pcount; pnum++) {
    sim->send_arrays[pnum].len = 0;
    sim->send_arrays[pnum].allocation = 0;

    sim->recv_arrays[pnum].len = 0;
    sim->recv_arrays[pnum].allocation = 0;
  }

  // datatypes
  MPI_Datatype type_vec3d;
  MPI_Type_contiguous(3, MPI_DOUBLE, &type_vec3d);
  MPI_Type_commit(&type_vec3d);

  // FRAGILE, any change to particle requires a change to this
  int block_lens[4] = {1, 1, 1, 1};
  MPI_Aint displacements[4] = {
      offsetof(Particle, position), offsetof(Particle, velocity),
      offsetof(Particle, acceleration), offsetof(Particle, mass)};
  MPI_Datatype types[4] = {type_vec3d, type_vec3d, type_vec3d, MPI_DOUBLE};
  MPI_Datatype type_particle;
  MPI_Type_create_struct(1, block_lens, displacements, types, &type_particle);
  MPI_Type_commit(&type_particle);
  sim->MPI_particle_type = type_particle;

  sim->use_mpi = true;
}

void MPI_teardown(Simulation *sim) {
  // TODO: free memory
  MPI_Finalize();
}
void MPI_sync_particles(Simulation *sim) {
  // throw out all non-owned particles
  // NOTE: a little hacky?
  sim->parray.len = sim->N_owned;

  // take turns sending size of payload
  int *recv_sizes = calloc(sim->MPI_pcount, sizeof(int));
  for (int pnum = 0; pnum < sim->MPI_pcount; pnum++) {
    if (pnum == sim->MPI_rank) {
      recv_sizes[pnum] = sim->N_owned;
    }
    MPI_Bcast(&recv_sizes[pnum], 1, MPI_INT, pnum, MPI_COMM_WORLD);
  }

  // allocate space for payload
  for (int pnum = 0; pnum < sim->MPI_pcount; pnum++) {
    if (pnum == sim->MPI_rank) {
      continue;
    }
    grow_array(&sim->recv_arrays[pnum], recv_sizes[pnum]);
  }

  // take turns recieving payload
  // TODO: should switch to non-blocking, major bottleneck

  //  setup non-blocking recieves
  MPI_Request *requests = calloc(sim->MPI_pcount, sizeof(MPI_Request));
  for (int pnum = 0; pnum < sim->MPI_pcount; pnum++) {
    if (pnum == sim->MPI_rank || recv_sizes[pnum] == 0) {
      continue;
    }
    MPI_Irecv(sim->recv_arrays[pnum].particles, recv_sizes[pnum],
              sim->MPI_particle_type, pnum, 0, MPI_COMM_WORLD, &requests[pnum]);
  }

  // send
  for (int pnum = 0; pnum < sim->MPI_pcount; pnum++) {
    if (pnum == sim->MPI_rank) {
      continue;
    }
    MPI_Send(sim->parray.particles, sim->parray.len, sim->MPI_particle_type,
             pnum, 0, MPI_COMM_WORLD);
  }

  // add to parray
  for (int pnum = 0; pnum < sim->MPI_pcount; pnum++) {
    if (pnum == sim->MPI_rank) {
      continue;
    }
    ParticleArray *n_parray = &sim->recv_arrays[pnum];
    n_parray->len = recv_sizes[pnum];
    for (int i = 0; i < n_parray->len; i++) {
      add_particle(sim, n_parray->particles[i], false);
    }
  }

  free(recv_sizes);
  free(requests);
}
