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
  // FRAGILE, any change to particle requires a change to this
  int p_block_lens[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Aint p_displacements[10] = {
      offsetof(Particle, pos_x), offsetof(Particle, pos_y),
      offsetof(Particle, pos_z), offsetof(Particle, vel_x),
      offsetof(Particle, vel_y), offsetof(Particle, vel_z),
      offsetof(Particle, acc_x), offsetof(Particle, acc_y),
      offsetof(Particle, acc_z), offsetof(Particle, mass)};
  printf("offsets: %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n",
         p_displacements[0], p_displacements[1], p_displacements[2],
         p_displacements[3], p_displacements[4], p_displacements[5],
         p_displacements[6], p_displacements[7], p_displacements[8],
         p_displacements[9]);
  MPI_Datatype p_types[10] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                              MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                              MPI_DOUBLE, MPI_DOUBLE};
  MPI_Datatype type_particle;
  MPI_Type_create_struct(1, p_block_lens, p_displacements, p_types,
                         &type_particle);
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

  printf("recv_sizes: %d %d %d %d\n", recv_sizes[0], recv_sizes[1],
         recv_sizes[2], recv_sizes[3]);

  // allocate space for payload
  for (int pnum = 0; pnum < sim->MPI_pcount; pnum++) {
    if (pnum == sim->MPI_rank) {
      continue;
    }
    grow_array(&sim->recv_arrays[pnum], recv_sizes[pnum]);
  }

  // take turns recieving payload
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
    MPI_Send(sim->parray.particles, sim->N_owned, sim->MPI_particle_type, pnum,
             0, MPI_COMM_WORLD);
  }

  // wait
  for (int pnum = 0; pnum < sim->MPI_pcount; pnum++) {
    if (pnum == sim->MPI_rank) {
      continue;
    }
    MPI_Status status;
    MPI_Wait(&requests[pnum], &status);
    int count;
    MPI_Get_count(&status, MPI_BYTE, &count);
    printf("recieved a count of %d\n", count);
  }

  // add to parray
  for (int pnum = 0; pnum < sim->MPI_pcount; pnum++) {
    if (pnum == sim->MPI_rank) {
      continue;
    }
    ParticleArray *n_parray = &sim->recv_arrays[pnum];
    n_parray->len = recv_sizes[pnum];
    for (int i = 0; i < n_parray->len; i++) {
      printf("(%.3f, %.3f, %.3f : %f)\n", n_parray->particles[i].vel_x,
             n_parray->particles[i].vel_y, n_parray->particles[i].vel_z,
             n_parray->particles[i].mass);
      add_particle(sim, n_parray->particles[i], false);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  free(recv_sizes);
  free(requests);
}
