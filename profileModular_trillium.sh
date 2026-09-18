#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --account=def-tobi
#SBATCH -t 00:25:00
#SBATCH -o profile_%j.out

# Real instruction/source-line profiling of the modular Boris mover
# (move_and_bc_boris, mod_particleMover.f90) on the ITER case, using
# `perf record` with DWARF call-graph unwinding.
#
# Context: source-level investigation (compiler vectorization reports,
# the B-field negligible-|B| skip ratio, and a full ParticleSet SoA->flat-
# array rewrite) has ruled out three plausible causes for modular's mover
# being ~2-4x slower than legacy's at ITER scale, without closing the gap
# any of those times. This script gets an actual instruction-level answer
# instead of a fourth source-level guess.
#
# Usage: sbatch profileModular_trillium.sh
#
# Output (read in this order):
#   profile_<jobid>.out              - this script's own log; check the
#                                       "perf permission check" section
#                                       first - if that failed, everything
#                                       below it did not run as intended.
#   perf_annotate_boris_<jobid>.txt  - per-source-line % of time inside
#                                       move_and_bc_boris specifically -
#                                       almost certainly the file to look
#                                       at first.
#   perf_report_flat_<jobid>.txt     - flat, self-time-sorted hotspot list
#                                       across the whole run (all phases,
#                                       not just the mover) for context.
#   perf_report_callgraph_<jobid>.txt - same data with caller call-graphs,
#                                        depth-limited to stay readable.
#
# EDIT BEFORE SUBMITTING: add whatever `module load` lines this cluster
# needs for the Intel oneAPI/MPI toolchain (mpiifx, mpirun) AND for `perf`
# if it's provided via a module here rather than always on PATH - unknown
# from here, matching the same TODO in runModular_trillium_calibrate.sh.

set -u
cd "$SLURM_SUBMIT_DIR"
JOBID="${SLURM_JOB_ID:-manual}"

echo "===== perf permission check ====="
# perf record needs perf_event_paranoid <= 1 (or CAP_PERFMON/root) to
# sample an unprivileged process's call stacks. Check this FIRST, before
# spending a build + a run on a profile that will come back empty or
# refuse to start - a HPC login/compute node commonly locks this down.
if [ -r /proc/sys/kernel/perf_event_paranoid ]; then
  paranoid=$(cat /proc/sys/kernel/perf_event_paranoid)
  echo "perf_event_paranoid = ${paranoid} (need <= 1 for call-graph sampling as a normal user)"
else
  echo "/proc/sys/kernel/perf_event_paranoid not readable - unknown restriction level"
fi
if ! command -v perf >/dev/null 2>&1; then
  echo "perf NOT FOUND on PATH - add a module load line for it above, or fall back to VTune"
  echo "  (manual VTune fallback: mpirun -np 8 vtune -collect hotspots -result-dir vtune_r%q{PMI_RANK} -- ./build_prof/run_min ,"
  echo "   then vtune -report hotspots -result-dir vtune_r0 -source-object function=move_and_bc_boris)"
  exit 1
fi
perf_ok=1
if ! timeout 10 perf record -o /tmp/perf_selftest_${JOBID}.data -- true >/tmp/perf_selftest_${JOBID}.log 2>&1; then
  perf_ok=0
fi
rm -f /tmp/perf_selftest_${JOBID}.data /tmp/perf_selftest_${JOBID}.log
if [ "$perf_ok" -eq 0 ]; then
  echo "perf record self-test FAILED - this node/account cannot sample call stacks as a normal user."
  echo "Ask cluster support to either raise perf_event_paranoid or grant CAP_PERFMON, or use the"
  echo "manual VTune fallback printed above instead."
  exit 1
fi
echo "perf record self-test OK - proceeding."
echo "======================================"
echo

echo "===== building modular (run_min) WITH debug info, separate build_prof/ dir ====="
# -g added on top of the normal Release flags (does not disable -O3/-xHost/
# etc, just keeps symbol table + line info) so perf/vtune can attribute
# samples to actual source lines instead of only bare function names.
# Separate build_prof/ directory so this never touches/invalidates the
# plain build/ directory runModular_trillium_calibrate.sh uses.
rm -rf build_prof
FC=mpiifx cmake -S . -B build_prof -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_FLAGS=-g
cmake --build build_prof --target run_min -j"$(nproc)"

if [ ! -x build_prof/run_min ]; then
  echo "BUILD FAILED - fix compile errors above before profiling."
  exit 1
fi

ulimit -s unlimited
mkdir -p DATA/DATA_2D
mkdir -p Output/Output_2D

# --- swap in the ITER case, with the two known fixes applied on top ---
# (same as runModular_trillium_calibrate.sh - see that script's comments)
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

echo "===== profiling: MPI=1 x OMP=16 (deliberately NOT the fastest split) ====="
echo "===== A local smoke test of this exact pipeline (8 threads, ~10s,     ==="
echo "===== -F 999) produced a 1.2GB perf.data file - --call-graph dwarf's  ==="
echo "===== per-sample stack dump dominates size, not just sample rate. The ==="
echo "===== full 192-thread/180s split from the calibrate sweep would scale ==="
echo "===== to hundreds of GB at that rate, so this run trades scale (fewer ==="
echo "===== threads, lower frequency, shorter cap) for a file that's       ==="
echo "===== actually possible to move/open - it doesn't change WHERE time  ==="
echo "===== goes inside a single call to move_and_bc_boris, which is the   ==="
echo "===== question this is answering (not which MPI/OMP split is fastest ==="
echo "===== - the sweeps already answered that separately)."
export OMP_NUM_THREADS=16 OMP_PROC_BIND=true OMP_PLACES=cores I_MPI_PIN_DOMAIN=omp

PERFDATA="perf_modular_1x16_${JOBID}.data"
# --call-graph dwarf: unwinds via the -g debug info instead of frame
# pointers, which -O3 may have omitted. perf record follows mpirun's
# forked/exec'd Hydra proxy and MPI rank children by default (no extra
# flag needed) - it does NOT need root/-a for that, only for true
# system-wide profiling of processes outside this one's descendant tree.
# -F 49: deliberately low frequency (default max is often 4000+) - see the
# file-size note above. timeout is a safety net; a SIGTERM here still lets
# perf flush a valid (if truncated) data file.
timeout 90 perf record -g --call-graph dwarf -F 49 -o "$PERFDATA" -- \
    mpirun -np 1 ./build_prof/run_min > "run_output_${JOBID}.log" 2>&1
rc=$?
echo "perf record + run exit code: $rc (124 = hit the 90s timeout - still fine, data up to that point is valid)"
echo "perf.data size: $(du -h "$PERFDATA" 2>/dev/null | cut -f1)"

pkill -9 -f "run_min" 2>/dev/null
pkill -9 -f "hydra_pmi_proxy" 2>/dev/null
pkill -9 -f "mpiexec.hydra" 2>/dev/null
sleep 2

if [ ! -s "$PERFDATA" ]; then
  echo "NO PERF DATA WRITTEN ($PERFDATA missing or empty) - check run_output_${JOBID}.log for a crash"
  echo "before perf even got samples (e.g. the run itself erroring immediately)."
  exit 1
fi

echo "===== generating reports from $PERFDATA ====="

# Flat, self-time-sorted hotspot list across the WHOLE run (poisson, mover,
# deposit, everything) - context for how big the mover's slice really is.
perf report --stdio --no-children -i "$PERFDATA" \
    > "perf_report_flat_${JOBID}.txt" 2>&1

# Same data with caller call-graphs, depth-limited so it stays readable
# rather than dumping every possible call chain.
perf report --stdio --call-graph=graph,3,caller -i "$PERFDATA" \
    > "perf_report_callgraph_${JOBID}.txt" 2>&1

# Full per-source-line annotation dump (can be large - full symbol table),
# then the move_and_bc_boris-specific slice pulled out into its own small,
# easy-to-paste-back file. mod_particleMover.f90 module procedures are
# name-mangled as mod_particlemover_mp_move_and_bc_boris_ (confirmed
# earlier from the ifx -qopt-report output), but grep on the unmangled
# substring "move_and_bc_boris" to not depend on that exact mangling.
perf annotate --stdio -i "$PERFDATA" > "perf_annotate_full_${JOBID}.txt" 2>&1

# Each hot symbol's disassembly in `perf annotate --stdio` output is its
# own block starting with a " Percent |...Disassembly of..." header line
# (NOT separated by blank lines - annotated source/asm lines routinely
# have a bare " :" with no percent, which looks blank-ish but isn't an
# empty line, so a blank-line-based extraction silently grabs the rest of
# the file). Buffer lines since the last such header; once a buffered
# block is found to contain move_and_bc_boris, print that whole block and
# stop at the next header.
awk '
  function flush() {
    if (active) { printf "%s", block; exit }
  }
  /^ Percent/ {
    flush()
    block = $0 "\n"
    active = 0
    next
  }
  { block = block $0 "\n" }
  /move_and_bc_boris/ { active = 1 }
  END { flush() }
' "perf_annotate_full_${JOBID}.txt" > "perf_annotate_boris_${JOBID}.txt"

if [ ! -s "perf_annotate_boris_${JOBID}.txt" ]; then
  echo "WARNING: perf_annotate_boris_${JOBID}.txt came out empty - the awk grab on"
  echo "'move_and_bc_boris' didn't find/bound the section as expected. The full dump"
  echo "(perf_annotate_full_${JOBID}.txt) still has it - grep that file by hand for"
  echo "'move_and_bc_boris' instead."
fi

echo "PROFILING DONE."
echo "Read perf_annotate_boris_${JOBID}.txt first (per-source-line %% inside the mover)."
echo "perf_report_flat_${JOBID}.txt for whole-run context, perf_report_callgraph_${JOBID}.txt for callers."
