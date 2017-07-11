#!/bin/bash

BIN=~/soft/phfmol_qmcpack/build/bin/qmcpack
ref="ref.xml"

for idet in {0..4};do
  cp $ref vmc.xml
  sed -i "s/det0.h5/det$idet.h5/" vmc.xml

  $BIN vmc.xml > out
  grep PASS out

  subdir=det$idet
  mkdir $subdir
  mv fftbox*.dat $subdir
  mv psig*.dat $subdir
done
