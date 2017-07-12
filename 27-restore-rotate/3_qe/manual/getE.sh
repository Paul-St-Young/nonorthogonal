#!/bin/bash
echo "#  ecut  energy walltime" > convergence.dat
for ecut in  20 40 80 160 320; do
  subdir=ecut$ecut
  E=`grep ! $subdir/scf.out | awk '{print $5}'`
  t=`grep PWSCF $subdir/scf.out | grep WALL | awk '{print $5}'`
  echo $ecut $E $t >> convergence.dat
done
