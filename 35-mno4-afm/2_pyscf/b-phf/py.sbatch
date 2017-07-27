#!/bin/bash
#SBATCH --job-name mno-dump
#SBATCH -N 1
#SBATCH --tasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH -t 00:30:00
#SBATCH -o phf.out
#SBATCH -e phf.err
#SBATCH -p pdebug

export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_avx.so:$MKLROOT/lib/intel64/libmkl_core.so:$MKLROOT/lib/intel64/libmkl_sequential.so:$LD_PRELOAD
export OMP_NUM_THREADS=36
srun -n 1 python run_phfmol.py
