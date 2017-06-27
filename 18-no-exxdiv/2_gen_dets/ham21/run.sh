#!/bin/bash

BIN=hmmt_phfmol.x
# to get this binary, modify hamilt.f90 and dump hmmt and ovmt in shutdown_hf
#  write in the same format as determinant to reuse parsing script
#  
#  write(*,*) 'Hamiltonian: 1'
#  write(*,*)  hmmt
#  write(*,*) 'Overlap:  1'
#  write(*,*)  ovmt

# link integral file
ln -s ../fcidump.dat
# copy over det.21
idet=21
cp ../phfrun.det.$idet phfrun.det
cp ../phfrun.dat.$idet phfrun.dat
# execute
$BIN > ham${idet}.out
