#!/bin/bash

BIN=phfmol.x

for k in {2..50}; do
  mv phfrun.out phfrun.out.$k
  cp phfrun.det phfrun.det.$k
  cp phfrun.dat phfrun.dat.$k
  mv fort.100 dets.$k 
  $BIN >& out
done
