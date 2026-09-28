#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=24
#SBATCH --account=def-tobi
#SBATCH -t 00:20:00
#SBATCH -o debuglegacy_%j.out

# Rebuilds legacy (3dphpic.exe) with bounds/uninit checking and full debug
# symbols, then reruns the ITER case just long enough to hit the SIGSEGV
# that job 2448619 hit reliably at it=1 (first collision episode) across
# all 4 splits tested. Goal: get an exact file:line instead of the
# "Unknown Unknown Unknown" unsymbolized traceback from the release build.
#
# -O0 (no -ipo/-xHost) is deliberate here, not an oversight: -CB (array
# bounds checking) and -check uninit are most reliable/precise without
# aggressive optimization reordering code around them - mirrors how
# CMakeLists.txt's own INTEL_DEBUG_FLAGS for the modular build drops
# -ipo too. This build will be MUCH slower than production - that's fine,
# it only needs to survive to the first collision call, not finish a
# real run.
#
# Usage: sbatch debugLegacyCrash_trillium.sh

set -u
cd "$SLURM_SUBMIT_DIR"

echo "===== building legacy (3dphpic.exe) WITH bounds/uninit checking ====="
# rm the submit-dir binary FIRST: the previous version of this script only
# checked "does 3dphpic.exe exist", which a failed build satisfies
# trivially if a stale binary from an earlier (different) build is still
# sitting here - that's exactly what happened last time (two real compile
# errors got masked because the script silently re-ran the OLD release
# binary and reported the same crash as if the debug build had run).
# Removing it first means a failed build leaves NO binary, so the check
# below is actually meaningful.
rm -f 3dphpic.exe
DEBUG_OPT="-O0 -g -traceback -CB -warn all,noexternal -check uninit -check pointers -check output_conversion -check format -fpe0 -qopenmp -auto"
( cd Src && make clean >/dev/null 2>&1; make OPT="$DEBUG_OPT" -j8 ) || \
( cd Src && make OPT="$DEBUG_OPT" )

if [ ! -x 3dphpic.exe ]; then
  echo "BUILD FAILED - fix compile errors above."
  exit 1
fi

# --- swap in the ITER case, same fixes the calibration script applies ---
cp input_dir/conditions.inp  input_dir/conditions.inp.bak
cp input_dir/geometry.inp    input_dir/geometry.inp.bak
cp input_dir/boundary.inp    input_dir/boundary.inp.bak
cp input_dir/ITER/conditions.inp input_dir/conditions.inp
cp input_dir/ITER/geometry.inp   input_dir/geometry.inp
cp input_dir/ITER/boundary.inp   input_dir/boundary.inp
sed -i 's/^-0.18 2 no 0\. 0\.2 0\.6 0\.8 !/-0.18 2 no 0. 0.2 0.6 0.8 0. !/' input_dir/conditions.inp
sed -i 's/^1000 ! frequency for saving data/5 ! frequency for saving data/' input_dir/conditions.inp

restore_inputs() {
  mv -f input_dir/conditions.inp.bak input_dir/conditions.inp
  mv -f input_dir/geometry.inp.bak   input_dir/geometry.inp
  mv -f input_dir/boundary.inp.bak   input_dir/boundary.inp
}
trap restore_inputs EXIT

ulimit -s unlimited
mkdir -p DATA/DATA_2D

echo "===== running debug legacy (8x24, bounds-checked) - expect it to be slow ====="
# OMP_STACKSIZE=2G, generous: -O0 debug builds use MORE stack per frame
# than the -O3 release build did, and OMP worker threads (where
# collision_OMP runs, inside !$OMP PARALLEL) get their own stack sized by
# this var - completely separate from ulimit -s above, which only covers
# the main thread. If this build no longer crashes, the original SIGSEGV
# was a stack-size problem, not a real out-of-bounds bug (-CB would have
# reported a specific bounds violation instead of a raw SIGSEGV if it
# were the latter).
env OMP_NUM_THREADS=24 OMP_STACKSIZE=2G OMP_PROC_BIND=true OMP_PLACES=cores I_MPI_PIN_DOMAIN=omp \
  timeout 900 mpirun -np 8 ./3dphpic.exe > "debuglegacy_run_${SLURM_JOB_ID}.log" 2>&1
rc=$?
echo "run exit code: $rc"
echo "===== tail of debuglegacy_run_${SLURM_JOB_ID}.log (the symbolized crash should be here) ====="
tail -80 "debuglegacy_run_${SLURM_JOB_ID}.log"
