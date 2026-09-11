#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --account=def-tobi
#SBATCH -t 00:30:00
#SBATCH -o calibrate_%j.out

# Trillium node: 2x AMD EPYC 9655 (Zen 5), 96 cores/socket, 192 cores/node.
# Purpose: this is NOT a production run. It answers one question before we
# commit to a production MPI x OMP split: how many real NUMA domains does
# this node actually expose, and which of a few candidate splits is fastest
# for the modular code specifically (its own per-rank OMP load balance was
# seen to degrade as OMP threads/rank grew past ~24-32 in our 2-socket,
# 32-core dev-box tests - untested at the thread counts a 96-core socket
# implies, so don't assume 2x96 is safe without checking).
#
# Usage: sbatch runModular_trillium_calibrate.sh
# Then read calibrate_<jobid>.out - it prints the NUMA topology first, then
# one short timed sample of the ITER case per candidate split below.

set -u
cd "$SLURM_SUBMIT_DIR"

echo "===== NUMA topology on this node ====="
lscpu | grep -E "^Socket|^NUMA|^Core|^Thread|^Model name"
echo
numactl --hardware 2>&1 || echo "(numactl not available - relying on lscpu above)"
echo "======================================"
echo

ulimit -s unlimited

# --- swap in the ITER case, with the two known fixes applied on top ---
cp input_dir/conditions.inp  input_dir/conditions.inp.bak
cp input_dir/geometry.inp    input_dir/geometry.inp.bak
cp input_dir/boundary.inp    input_dir/boundary.inp.bak
cp input_dir/ITER/conditions.inp input_dir/conditions.inp
cp input_dir/ITER/geometry.inp   input_dir/geometry.inp
cp input_dir/ITER/boundary.inp   input_dir/boundary.inp
# missing trailing E0_RF field in the example file - without this the input
# parser desyncs and read_input aborts with "Input file was not read correctly!"
sed -i 's/^-0.18 2 no 0\. 0\.2 0\.6 0\.8 !/-0.18 2 no 0. 0.2 0.6 0.8 0. !/' input_dir/conditions.inp
# print diagnostics every 5 steps instead of 1000, so a short calibration
# run actually produces a few timing blocks to look at
sed -i 's/^1000 ! frequency for saving data/5 ! frequency for saving data/' input_dir/conditions.inp

restore_inputs() {
  mv input_dir/conditions.inp.bak input_dir/conditions.inp
  mv input_dir/geometry.inp.bak   input_dir/geometry.inp
  mv input_dir/boundary.inp.bak   input_dir/boundary.inp
}
trap restore_inputs EXIT

run_split() {
  local mpi_ranks="$1" omp_threads="$2"
  echo "--- MPI=${mpi_ranks} x OMP=${omp_threads} (= $((mpi_ranks*omp_threads)) cores) ---"
  env OMP_NUM_THREADS="$omp_threads" OMP_PROC_BIND=true OMP_PLACES=cores \
      I_MPI_PIN_DOMAIN=omp \
      timeout 90 mpirun -np "$mpi_ranks" ./build/run_min \
      2>&1 | grep -E "TIME STEP|total      \(ms\)|mover      \(ms\)|poisson    \(ms\)|max/avg"
  echo
}

# Candidates to compare - adjust once the topology above tells you the real
# NUMA domain count. These are reasonable starting guesses for a 2-socket,
# 96-cores/socket node under a few common NPS (NUMA-per-socket) BIOS settings:
run_split 2  96   # 1 rank per socket, NPS1 assumption - baseline, matches our validated 2x16 pattern scaled up
run_split 6  32   # 3 ranks per socket - hedges against OMP imbalance growth per rank
run_split 8  24   # matches NPS4 (4 NUMA domains/socket) if that's how this node is configured
run_split 4  48   # matches NPS2 (2 NUMA domains/socket)

echo "SWEEP DONE - compare the 'total (ms)' steady-state values above."
