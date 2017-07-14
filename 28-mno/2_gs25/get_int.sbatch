#!/bin/bash
#SBATCH --job-name dump
#SBATCH -N 1
#SBATCH --tasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH -t 00:10:00
#SBATCH -o int.out
#SBATCH -e int.err
#SBATCH -p pdebug

# use MKL libraries, this is critical for speed! (2min MKL -> 30min numpy)
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_avx.so:$MKLROOT/lib/intel64/libmkl_core.so:$MKLROOT/lib/intel64/libmkl_sequential.so:$LD_PRELOAD
export OMP_NUM_THREADS=36
srun -n 1 python run_int.py
