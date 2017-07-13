#!/bin/bash
#SBATCH --job-name pyscf
#SBATCH -N 1
#SBATCH --tasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH -t 00:10:00
#SBATCH -o default.out
#SBATCH -e default.err
#SBATCH -p pbatch

export OMP_NUM_THREADS=36
srun -n 1 python main.py
