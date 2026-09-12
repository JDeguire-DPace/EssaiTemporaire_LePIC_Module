#!/bin/bash
# Local-box equivalent of runModular_trillium_calibrate.sh - NOT a Slurm job,
# just run it directly: ./runModular_local_calibrate.sh
#
# This dev box: 2x Intel Xeon Gold 6346 @ 3.10GHz, 16 physical cores/socket,
# 2 threads/core (HT on) -> 64 logical CPUs, but only 32 DISTINCT physical
# cores. lscpu/cgroup confirm 2 NUMA nodes:
#   node0: logical CPUs 0-15  (physical) + 32-47 (their HT siblings)
#   node1: logical CPUs 16-31 (physical) + 48-63 (their HT siblings)
# i.e. logical IDs 0-31 ARE the 32 physical cores, one per id, no overlap;
# 32-63 are just the hyperthread half of 0-31. So every candidate split below
# totals to 32 -- matching how the Trillium script always totalled 192 -- and
# OMP_PLACES=cores keeps each OMP thread on its own physical core (avoiding
# HT siblings) automatically, no numactl/taskset pinning needed for that.
#
# Purpose: same two questions as the cluster version, sized for 2 domains of
# 16 cores instead of 8 domains of 24:
#   1. Confirm this box really is a flat 2-domain (1 domain/socket) layout.
#   2. Which MPI x OMP split at 32 physical cores is fastest, modular vs
#      legacy - plus a bonus HT-enabled data point (2x32, using all 64
#      logical CPUs) since this box, unlike Trillium's, actually has SMT.
#
# NOTE: mpiifx/mpirun (Intel oneAPI 2021.14) are already on PATH here - no
# module system needed locally, unlike the cluster version.

set -u
cd "$(dirname "$(readlink -f "$0")")"

echo "===== NUMA topology on this box ====="
lscpu | grep -E "^Socket|^NUMA|^Core|^Thread|^Model name|^CPU\(s\):"
echo
numactl --hardware 2>&1 || echo "(numactl not available - relying on lscpu above; physical-core IDs are 0-31, see header comment)"
echo "======================================"
echo

# --- build both executables fresh, on this box/architecture ---
# Same reasoning as the cluster script: don't trust a build/ that might have
# been configured differently earlier - rebuild clean so -xHost is correct
# for THIS run.
echo "===== building modular (run_min) ====="
rm -rf build
FC=mpiifx cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target run_min -j"$(nproc --all)"

echo "===== building legacy (3dphpic.exe) ====="
# Same known Makefile race as on the cluster if -j is too high for ~20 small
# source files - use a modest -j and retry once serially if it still fails.
( cd Src && make clean >/dev/null 2>&1; make -j8 ) || \
( cd Src && make )

if [ ! -x build/run_min ] || [ ! -x 3dphpic.exe ]; then
  echo "BUILD FAILED - fix compile errors above before calibrating."
  exit 1
fi

ulimit -s unlimited
mkdir -p DATA/DATA_2D

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
  mv -f input_dir/conditions.inp.bak input_dir/conditions.inp
  mv -f input_dir/geometry.inp.bak   input_dir/geometry.inp
  mv -f input_dir/boundary.inp.bak   input_dir/boundary.inp
}
trap restore_inputs EXIT

run_split() {
  local label="$1" exe="$2" mpi_ranks="$3" omp_threads="$4" omp_places="${5:-cores}"
  echo "--- ${label}: MPI=${mpi_ranks} x OMP=${omp_threads} (= $((mpi_ranks*omp_threads)) cores, OMP_PLACES=${omp_places}) ---"
  env OMP_NUM_THREADS="$omp_threads" OMP_PROC_BIND=true OMP_PLACES="$omp_places" \
      I_MPI_PIN_DOMAIN=omp \
      timeout 90 mpirun -np "$mpi_ranks" "$exe" \
      2>&1 | grep -E "TIME STEP|total      \(ms\)|mover      \(ms\)|poisson    \(ms\)|max/avg|^ it=|^ <t>"
  echo
  # Same zombie-rank safety net as the cluster script.
  pkill -9 -f "run_min" 2>/dev/null
  pkill -9 -f "3dphpic.exe" 2>/dev/null
  pkill -9 -f "hydra_pmi_proxy" 2>/dev/null
  pkill -9 -f "mpiexec.hydra" 2>/dev/null
  sleep 2
}

# Confirmed from this box's own lscpu/cgroup output: 2 NUMA domains of 16
# physical cores each (1 domain/socket - unlike Trillium's 4 domains/socket).
# 2x16 is the one candidate that lines up exactly with real hardware; 1x32
# spans both domains (the whole box) in a single rank; 4x8 and 8x4 test
# going finer than a domain, by increasing amounts.
for split in "2 16" "1 32" "4 8" "8 4"; do
  read -r ranks threads <<< "$split"
  run_split "modular" "./build/run_min" "$ranks" "$threads"
  run_split "legacy"  "./3dphpic.exe"   "$ranks" "$threads"
done

# Bonus: this box actually has hyperthreading (unlike the Trillium node),
# so also check whether using both HT siblings per core helps or hurts.
# OMP_PLACES=threads here (not "cores") so each of the 32 OMP threads/rank
# gets its own hardware thread instead of OMP_PLACES=cores failing to find
# 32 distinct physical cores inside a 16-physical-core NUMA domain.
run_split "modular (HT)" "./build/run_min" 2 32 threads
run_split "legacy (HT)"  "./3dphpic.exe"   2 32 threads

echo "SWEEP DONE - compare the 'total (ms)' steady-state values above (skip any block where MC/sorting just fired - those spike and aren't representative of a typical step)."
