#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=48
#SBATCH --account=def-tobi
#SBATCH -t 22:45:01

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
export OMP_PLACES=cores

echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

mpirun -np 2 ./build/run_min > dump.modular
