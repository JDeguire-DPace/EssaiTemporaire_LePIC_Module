#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=48
#SBATCH --account=def-tobi
#SBATCH -t 1:00:01

export OMP_NUM_THREADS=48
export OMP_PROC_BINF=true

mpirun -np 4 ./run_min > run.dump
