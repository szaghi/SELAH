#!/usr/bin/env bash
# Run SELAH's verification programs in report mode and regenerate the documentation's data, tables, figures and casts.
#
# Usage: scripts/run_vv.sh [--mode MODE] [--data DIR] [--no-build] [--no-figures] [--only NAME ...]
#   --mode MODE    fobis mode that builds the programs (default: tests-gnu)
#   --data DIR     where the report data goes (default: exe/vv-data; build products, not committed)
#   --no-build     use the executables already in exe/
#   --no-figures   run the programs only
#   --only NAME    run only these programs; repeatable
#
# Every verification program asserts its bounds in both modes; with SELAH_VV_DIR set it also writes its data there.
# The figures (light and dark twins, docs/scripts/figure_theme.py) and the animations go to docs/public/verification/,
# the result tables to docs/verification/_tables/ (included by the case pages). Needs a Python with numpy and
# matplotlib (SELAH_VV_PYTHON, default exe/.vv-venv/bin/python, else python3) and ffmpeg for the animations.
set -euo pipefail
cd "$(dirname "$0")/.."
MODE=tests-gnu
DATA=exe/vv-data
BUILD=1
FIGURES=1
ONLY=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --data) DATA="$2"; shift 2 ;;
    --no-build) BUILD=0; shift ;;
    --no-figures) FIGURES=0; shift ;;
    --only) ONLY+=("$2"); shift 2 ;;
    *) printf 'unknown argument: %s\n' "$1" >&2; exit 2 ;;
  esac
done
PY="${SELAH_VV_PYTHON:-exe/.vv-venv/bin/python}"
[[ -x $PY ]] || PY=python3
[[ $BUILD == 1 ]] && fobis build --mode "$MODE"
PROGRAMS=(selah_test_schemes selah_test_splittings selah_test_stepper selah_vv_stability selah_vv_ssp_tvd
          selah_vv_splitting_energy selah_vv_forms)
[[ ${#ONLY[@]} -gt 0 ]] && PROGRAMS=("${ONLY[@]}")
mkdir -p "$DATA"
DATA_ABS="$(cd "$DATA" && pwd)"
status=0
for p in "${PROGRAMS[@]}"; do
  printf '%-28s ' "$p"
  start=$(date +%s)
  if SELAH_VV_DIR="$DATA_ABS" "./exe/$p" > "$DATA_ABS/$p.log" 2>&1; then
    printf 'ok    %5ss\n' "$(( $(date +%s) - start ))"
  else
    printf 'FAIL  %5ss  (see %s)\n' "$(( $(date +%s) - start ))" "$DATA/$p.log"; status=1
  fi
done
if [[ $FIGURES == 1 ]]; then
  "$PY" docs/scripts/make_doc_figures.py --data "$DATA_ABS" --out docs/public/verification --tables docs/verification/_tables
  "$PY" docs/scripts/make_doc_animations.py --data "$DATA_ABS" --out docs/public/verification/casts
fi
exit $status
