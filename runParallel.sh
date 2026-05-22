#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --account=def-tobi
#SBATCH -t 12:25:01

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
export OMP_PLACES=cores

echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

srun ./3dphpic.exe > run.dump
