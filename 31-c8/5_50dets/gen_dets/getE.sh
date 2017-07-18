#!/bin/bash

ndet=50
for k in `seq 0 $ndet`;do
  fname="phfrun.out.$k"
  E=`grep 'f     =' $fname | tail -n 1 | awk '{print $3-6.673036300700005}'`
  echo $k $E >> ewald.dat #ke.dat
done

#awk '{print $2-6.673036300700005}' ke.dat
