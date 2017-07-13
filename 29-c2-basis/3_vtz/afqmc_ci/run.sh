#!/bin/bash
#SBATCH -N 1
#SBATCH -t 10

export OMP_NUM_THREADS=1
BIN=/g/g91/yang41/soft/master_qmcpack/afqmc_build/bin/qmcpack
ndet=`ls ../gen_dets | grep mat | wc -l`

for idet in `seq 45 $ndet`; do
  cp ref.xml c2.xml
  sed -i "s/mydets/..\/gen_dets\/dets.$idet/" c2.xml
  $BIN c2.xml > out$idet
done
