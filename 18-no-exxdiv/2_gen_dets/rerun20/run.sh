#!/bin/bash

BIN=~/soft/phfmol/phfmol.x
# print mdvec before the second call print_darr_ao in phfscf.f90

# link integral file
ln -s ../fcidump.dat
idet=20
cp ../phfrun.det.$idet phfrun.det
cp ../phfrun.dat.$idet phfrun.dat
# execute
$BIN > ham${idet}.out
