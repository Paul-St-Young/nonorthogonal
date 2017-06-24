#!/bin/bash

BIN=phfmol.x

#663379
integral_file=../1_dump_integrals/fcidump.dat
inp0_ref=phfrun0.inp
inp1_ref=phfrun1.inp
ndet=50

# run first calculation to reproduce HF determinant
#  fort.100, phfrun.out, phfrun.dat and phfrun.det will be generated
ln -s $integral_file fcidump.dat
cp $inp0_ref phfrun.inp
$BIN >& err
mv phfrun.out phfrun.out.0
mv fort.100 dets.0

# run more calculations to generate a multi-determinant expansion
for k in `seq 1 $ndet`; do
  cp phfrun.det phfrun.det.$k
  cp phfrun.dat phfrun.dat.$k
  cp $inp1_ref phfrun.inp
  $BIN >& err
  mv phfrun.out phfrun.out.$k
  mv fort.100 dets.$k 
done
