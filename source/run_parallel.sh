#!/bin/bash

STEPS=10000
echo "Running experiments for $STEPS steps "
for((p=4; p <= 256; p*=2))
do
  for((n=256; n <= 4096; n=n*2))
  do
    echo " " 
    echo "processes: $p, particles: $n" 
    for ((j=1; j <= 5; j++)) 
    do
      echo "   run $j"
      time mpirun -np $p -hostfile hostfilehuge ./benchmark_mpi $n $STEPS
    done
  done
done

