#!/bin/bash
#SBATCH --job-name qe
#SBATCH -N 1
#SBATCH --tasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH -t 00:05:00
#SBATCH -o scf.out
#SBATCH -e scf.err
#SBATCH -p pdebug

export OMP_NUM_THREADS=1

for ecut in 20 40 80 160 320; do
  subdir=ecut$ecut
  cp -r ref $subdir
  cd $subdir
  sed -i "s/myecut/$ecut/" scf.in
  pw.x < scf.in > scf.out
  cd ..
done
