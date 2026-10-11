#!/usr/bin/env bash
#
# Reproduce every result in the README.
#
#   ./reproduce.sh
#
# Runtime: well under a minute on a 2-core laptop. CPU only, no GPU, no
# network access, no external data. Every figure in the README comes from
# these four scripts and nothing else.
#
set -euo pipefail

cd "$(dirname "$0")"

PY="${PYTHON:-python3}"

echo "=============================================================="
echo " yamanaka-replay -- full reproduction"
echo "=============================================================="
"$PY" --version
"$PY" -c "import numpy; print('numpy', numpy.__version__)"
echo

mkdir -p results

START=$SECONDS

"$PY" experiments/01_main_comparison.py
"$PY" experiments/02_prior_robustness.py
"$PY" experiments/03_misspecification.py
"$PY" experiments/04_trace_single_run.py

echo "=============================================================="
echo " done in $((SECONDS - START))s -- CSV outputs in results/"
echo "=============================================================="
