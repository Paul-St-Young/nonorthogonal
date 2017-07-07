#!/bin/bash

BIN=~/soft/git_phfmol/phfmol.x
$BIN

cat fort.12 >> phfrun.out
# ditch fort.100
mv phfrun.out det01.out
