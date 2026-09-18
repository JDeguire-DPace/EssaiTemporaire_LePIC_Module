#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --account=def-tobi
#SBATCH -t 00:45:00
#SBATCH -o calibrate_%j.out

# Trillium node: 2x AMD EPYC 9655 (Zen 5), 96 cores/socket, 192 cores/node.
# Purpose: this is NOT a production run. It answers two questions before we
# commit to a production MPI x OMP split:
#   1. How many real NUMA domains does this node actually expose?
#   2. Which of a few candidate splits is fastest, for BOTH codes - modular's
#      own per-rank OMP load balance was seen to degrade as OMP threads/rank
#      grew past ~24-32 in our 2-socket, 32-core dev-box tests, untested at
#      the thread counts a 96-core socket implies, so don't assume 2x96 is
#      safe without checking.
#
# Usage: sbatch runModular_trillium_calibrate.sh
#    or: sbatch --export=ALL,SPLIT_TIMEOUT=300 runModular_trillium_calibrate.sh
#        (per-split time budget, default 150s - this site's sbatch defaults
#        to --export=NONE, i.e. it does NOT inherit your shell env, so
#        `SPLIT_TIMEOUT=300 sbatch ...` silently has no effect here: you
#        must pass --export=ALL,SPLIT_TIMEOUT=... explicitly. If you raise
#        SPLIT_TIMEOUT a lot, also raise the #SBATCH -t above accordingly -
#        4 splits x 2 codes x SPLIT_TIMEOUT is the worst-case sweep time.)
# Then read calibrate_<jobid>.out - topology first, then one short timed
# sample of the ITER case per candidate split, for legacy and modular both.
# Each split also writes its full, unfiltered stdout+stderr to
# fulllog_<jobid>_<modular|legacy>_<ranks>x<threads>.log in the submit dir -
# calibrate_<jobid>.out itself only ever shows a grep'd summary plus a clear
# completed/CRASHED/TIMED OUT verdict per split, so check the matching
# fulllog_* file for the full picture (especially on CRASHED).
#
# EDIT BEFORE SUBMITTING: add whatever `module load` lines this cluster
# needs for the Intel oneAPI/MPI toolchain (mpiifx, mpirun) - unknown from
# here. Match whatever module gave you the compiler in the build step below.

set -u
cd "$SLURM_SUBMIT_DIR"

echo "===== NUMA topology on this node ====="
lscpu | grep -E "^Socket|^NUMA|^Core|^Thread|^Model name"
echo
numactl --hardware 2>&1 || echo "(numactl not available - relying on lscpu above)"
echo "======================================"
echo

# --- build both executables fresh, on this node/architecture ---
# Do NOT reuse a build/ or Src/*.o copied from another machine: the
# Makefile/CMake build uses -xHost, which bakes in the BUILD machine's
# detected ISA. A binary built on an Intel dev box may not run correctly
# (or may run un-vectorized) on this Zen 5 node - rebuild here.
echo "===== building modular (run_min) ====="
rm -rf build
# The "Fortran compiler basename"/"Using MPI wrapper" lines printed by this
# configure step are the check that matters - confirm they say mpiifx/TRUE
# in calibrate_<jobid>.out before trusting anything downstream.
FC=mpiifx cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target run_min -j"$(nproc)"

echo "===== building legacy (3dphpic.exe) ====="
# -j"$(nproc)" (192-way) on ~20 small source files was oversubscribed and
# triggered a known Makefile race (link step firing before every .o exists,
# "cp: cannot stat '3dphpic.exe'") - seen on the dev box too. Use a modest
# -j and retry once serially if it still fails.
( cd Src && make clean >/dev/null 2>&1; make -j8 ) || \
( cd Src && make )

if [ ! -x build/run_min ] || [ ! -x 3dphpic.exe ]; then
  echo "BUILD FAILED - fix compile errors above before calibrating."
  exit 1
fi

ulimit -s unlimited
# DATA/DATA_2D (no leading dot/Output) is legacy's own output convention;
# Output/Output_2D is the modular code's (see mod_simulation.f90,
# mod_generateBoundary.f90, mod_injection.f90) - both are needed since we
# run both executables. Output/Output_2D is normally git-tracked via a
# placeholder file so a real `git clone` already has it, but mkdir -p here
# too in case this scratch checkout didn't preserve that.
mkdir -p DATA/DATA_2D
mkdir -p Output/Output_2D

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

# How long each split gets before being killed. Override with
# SPLIT_TIMEOUT=<seconds> sbatch ...  if the ITER case turns out to need
# longer than this for a single step (large-scale cases seen so far run
# into the hundreds of millions of particles - step 1 in particular also
# eats particle loading + an initial full Poisson solve, so it can be much
# slower than the steady-state steps that follow).
: "${SPLIT_TIMEOUT:=150}"

run_split() {
  local label="$1" exe="$2" mpi_ranks="$3" omp_threads="$4"
  local logfile="fulllog_${SLURM_JOB_ID}_${label}_${mpi_ranks}x${omp_threads}.log"
  echo "--- ${label}: MPI=${mpi_ranks} x OMP=${omp_threads} (= $((mpi_ranks*omp_threads)) cores) - full output: ${logfile} ---"
  # Previous version of this script piped stdout+stderr straight through a
  # grep filter and discarded mpirun's real exit code - a crash (MPI_Abort,
  # segfault, OOM kill) produced a message that didn't match the filter, so
  # it silently vanished and looked identical to "still computing". Always
  # capture everything to a file first, THEN filter/report from that, so a
  # real error is never lost.
  env OMP_NUM_THREADS="$omp_threads" OMP_PROC_BIND=true OMP_PLACES=cores \
      I_MPI_PIN_DOMAIN=omp \
      timeout "$SPLIT_TIMEOUT" mpirun -np "$mpi_ranks" "$exe" \
      > "$logfile" 2>&1
  local rc=$?
  if [ $rc -eq 124 ]; then
    echo "  >>> TIMED OUT after ${SPLIT_TIMEOUT}s (was still running - may just be slow at this split, not broken; raise SPLIT_TIMEOUT to check) <<<"
  elif [ $rc -ne 0 ]; then
    echo "  >>> CRASHED, exit=${rc} - last lines of ${logfile}: <<<"
    tail -20 "$logfile" | sed 's/^/  | /'
  else
    echo "  >>> completed normally, exit=0 <<<"
  fi
  grep -E "TIME STEP|total      \(ms\)|mover      \(ms\)|poisson    \(ms\)|max/avg|^ it=|^ <t>|boris used" "$logfile"
  echo
  # Safety net: earlier testing found that killing/timing-out the mpirun
  # wrapper does not reliably reap the actual worker ranks underneath it
  # (Hydra's process tree), leaving zombies that silently contend with the
  # NEXT split's measurement. Force-clear everything between runs.
  pkill -9 -f "run_min" 2>/dev/null
  pkill -9 -f "3dphpic.exe" 2>/dev/null
  pkill -9 -f "hydra_pmi_proxy" 2>/dev/null
  pkill -9 -f "mpiexec.hydra" 2>/dev/null
  sleep 2
}

# Confirmed from this node's own numactl output: 8 NUMA domains of 24 cores
# each (4/socket, distance 12 within a socket's domains, 32 across sockets).
# 8x24 is the one candidate that lines up exactly with real hardware; 4x48
# still stays within one socket (2 domains/rank); 2x96 spans a whole socket
# (all 4 of its domains) per rank; 16x12 tests going finer than a domain.
for split in "8 24" "4 48" "2 96" "16 12"; do
  read -r ranks threads <<< "$split"
  run_split "modular" "./build/run_min" "$ranks" "$threads"
  run_split "legacy"  "./3dphpic.exe"   "$ranks" "$threads"
done

echo "SWEEP DONE - compare the 'total (ms)' steady-state values above (skip any block where MC/sorting just fired - those spike and aren't representative of a typical step)."
