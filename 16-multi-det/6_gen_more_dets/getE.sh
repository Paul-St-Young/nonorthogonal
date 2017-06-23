#!/bin/bash

ndet=`ls dets* | wc -l`

for k in `seq 0 $(($ndet-1))`; do

  energy=`grep 'f     =' phfrun.out.$k | tail -n 1 | awk '{print $3}'`
  echo "$k  $energy" >> ke.dat

done
