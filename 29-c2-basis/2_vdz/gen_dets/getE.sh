#!/bin/bash

ndet=50
for k in `seq 0 $ndet`;do
  fname="phfrun.out.$k"
  E=`grep 'f     =' $fname | tail -n 1 | awk '{print $3}'`
  echo $k $E >> ke.dat
done

awk '{print $2+7.562112440510}' ke.dat
