#!/bin/bash

STEPS=10000
echo "Running experiments for $STEPS steps "
for((n=256; n <= 8192; n=n*2))
do
  echo " " 
  echo "particles: $n" 
  for ((j=1; j <= 5; j++)) 
  do
    echo "   run $j"
    time ./benchmark_tree $n $STEPS
  done
done

