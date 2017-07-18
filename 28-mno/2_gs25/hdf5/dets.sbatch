#!/bin/bash
#SBATCH --job-name gen_dets
#SBATCH -N 1
#SBATCH -t 30
#SBATCH -p pdebug

BIN=/g/g91/yang41/notes/mmorale/hdf5_phfmol/phfmol.x
export OMP_NUM_THREADS=36

cp ../phfrun.det .
$BIN
