# CS87-Final-Project


This is where we keep our code. 

Maybe text documents.

## Results?
Here's the first results of running the particle-particle algorithm with and
without MPI. Promising stuff.

``` sh
$ time ./benchmark 1000 10000

real	0m28.828s
user	0m28.822s
sys	0m0.005s
$ time mpirun -np 4 ./benchmark_mpi 1000 10000

real	0m8.655s
user	0m33.896s
sys	0m0.110s
```
