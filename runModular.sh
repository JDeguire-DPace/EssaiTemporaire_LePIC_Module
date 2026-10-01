#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=24
#SBATCH --account=def-tobi
#SBATCH -t 23:45:01

# 8 ranks x 24 OMP threads/rank (= 192 cores, one rank per NUMA domain on
# Trillium's 8x24 topology per runModular_trillium_calibrate.sh's numactl
# output), NOT a single big-thread-count rank.
#
# Why: measured locally (2-socket/32-physical-core dev box, same ITER case)
# that modular's own per-rank OMP load balance is fine up to its physical
# core count but gets volatile and erratic past it (mover_push max/avg rose
# from ~1.0 at 8-32 threads/rank to spikes of 1.67 at 64 - three thread
# counts past this box's 32 physical cores), and separately that modular's
# per-thread scaling efficiency falls behind legacy's well before that even
# at "safe" thread counts (legacy: ~97% scaling efficiency 8->64 threads;
# modular: ~50%, plateauing hard past 32). Legacy showed neither effect
# (flat ~1.04-1.10 imbalance at every thread count tested, up to 64).
# Splitting the node into more, smaller-thread-count ranks avoids both:
# every rank's own OMP team stays comfortably inside the region where
# modular scales well, and MPI (not OpenMP) does the work of filling the
# rest of the node.
#
# This 8x24 split is carried over from the calibration sweep's NUMA-aligned
# candidate, NOT independently re-verified at ITER scale with multiple MPI
# ranks (today's local testing only exercised single-rank/many-OMP-thread
# configs) - rerun runModular_trillium_calibrate.sh on Trillium once
# available to confirm/tune this against the other candidate splits before
# treating it as final.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
# OMP_PLACES=threads, not cores: an overnight affinity sweep (OMP=48,
# oversubscribed on the dev box) measured threads giving the best imbalance
# ratio (1.091) vs cores (1.113) and sockets (1.84-2.06, much worse - avoid
# entirely). Difference is modest at cores vs threads but free to take.
export OMP_PLACES=threads
export I_MPI_PIN_DOMAIN=omp

echo "OMP_NUM_THREADS=$OMP_NUM_THREADS, MPI ranks=$SLURM_NTASKS"

mpirun -np "$SLURM_NTASKS" ./build/run_min > dump.modular
