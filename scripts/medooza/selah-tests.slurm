#!/bin/bash
# SELAH on medooza (CNR-IAC): build a test mode on a compute node and run every test program on an A30.
#
# The source tree is built where it lies (default /home/szaghi/selah, the build area) and the programs run from
# their own run directory on /ifs (default /ifs/hpc/home/szaghi/selah-tests/<job id>): nothing is written next to
# the sources but the build products. Submit from the login node (only sbatch runs there):
#
#   sbatch scripts/medooza/selah-tests.slurm                                  # tests-nvf-oac, varset iac
#   sbatch --export=ALL,SELAH_SRC=/home/szaghi/selah-m3 scripts/medooza/selah-tests.slurm
#   sbatch --export=ALL,SELAH_MODE=bench-nvf-oac,SELAH_RUN_GLOB='selah_bench_*' scripts/medooza/selah-tests.slurm
#
# The node drivers (545 on prod-gn, 575 on dev-gn) cannot run CUDA 13: varset iac builds for cc80 with CUDA 12.9.
#SBATCH -A iacroma
#SBATCH -p dev-gn
#SBATCH -J selah-tests
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH -t 00:50:00
#SBATCH -o /ifs/hpc/home/szaghi/selah-tests/selah-%j.out
#SBATCH -e /ifs/hpc/home/szaghi/selah-tests/selah-%j.out
set -uo pipefail

SRC=${SELAH_SRC:-/home/szaghi/selah}
MODE=${SELAH_MODE:-tests-nvf-oac}
VARSET=${SELAH_VARSET:-iac}
RUN_GLOB=${SELAH_RUN_GLOB:-selah_test_*}
RUN=${SELAH_RUN:-/ifs/hpc/home/szaghi/selah-tests/$SLURM_JOB_ID}

source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh || exit 1
export PATH=$HOME/.local/bin:$PATH
mkdir -p "$RUN"
{
   echo "job $SLURM_JOB_ID on $(hostname), $(date -Is)"
   echo "source $SRC, mode $MODE, varset $VARSET"
   nvfortran --version | sed -n 2p
   nvidia-smi --query-gpu=name,driver_version,clocks.max.sm,memory.total --format=csv,noheader
} | tee "$RUN/env.txt"

cd "$SRC" || exit 1
if ! fobis build --mode "$MODE" --varset "$VARSET" > "$RUN/build.log" 2>&1; then
   echo "BUILD FAILED: see $RUN/build.log" | tee "$RUN/summary.txt"
   exit 1
fi
# the build directory of the mode (fobos: tests-* build into exe/, bench-* into exe-bench/)
if [[ $MODE == bench-* ]]; then exe_dir=exe-bench; else exe_dir=exe; fi

cd "$RUN" || exit 1
for t in "$SRC/$exe_dir"/$RUN_GLOB; do
   [ -x "$t" ] && [ -f "$t" ] || continue
   name=$(basename "$t")
   if [[ $name == *_xfail_* ]]; then # expected failure: must exit non-zero (scripts/run_tests.sh convention)
      if "$t" > "$name.log" 2>&1; then echo "FAIL  $name (XPASS: exited 0)"; else echo "PASS  $name (XFAIL)"; fi
   elif "$t" > "$name.log" 2>&1; then
      echo "PASS  $name"
   else
      echo "FAIL  $name (exit $?)"
   fi
done | tee "$RUN/summary.txt"
if grep -q '^FAIL' "$RUN/summary.txt" || ! grep -q '^PASS' "$RUN/summary.txt"; then
   echo "some program failed, or none ran" | tee -a "$RUN/summary.txt"
   exit 1
fi
echo "all passed" | tee -a "$RUN/summary.txt"
