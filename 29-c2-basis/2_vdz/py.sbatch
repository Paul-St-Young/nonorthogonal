#!/bin/bash
#SBATCH --job-name pyscf
#SBATCH -N 1
#SBATCH --tasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH -t 00:10:00
#SBATCH -o default.out
#SBATCH -e default.err
#SBATCH -p pbatch

# use MKL libraries, this is critical for speed! (2min MKL -> 30min numpy)
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_avx.so:$MKLROOT/lib/intel64/libmkl_core.so:$MKLROOT/lib/intel64/libmkl_sequential.so:$LD_PRELOAD
export OMP_NUM_THREADS=36
srun -n 1 python main.py
