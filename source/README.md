# CS87-Final-Project


This is where we keep our code. 

Maybe text documents.

## Results?
Here are all the experiments we've run thus far, looks not bad at all.

``` sh
$ time ./benchmark 1000 10000

real	0m28.828s
user	0m28.822s
sys	0m0.005s
$ time mpirun -np 4 ./benchmark_mpi 1000 10000

real	0m8.655s
user	0m33.896s
sys	0m0.110s
$ time ./benchmark_tree 1000 10000

real	0m14.812s
user	0m14.487s
sys	0m0.324s
```

