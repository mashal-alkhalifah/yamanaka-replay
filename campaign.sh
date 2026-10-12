#!/usr/bin/env bash
#
# Run the large-scale campaign: experiments 05-10, 6,320 independent replicates.
#
#   ./campaign.sh
#
# This is NOT the reproduction path -- use ./reproduce.sh for that, which
# regenerates every README figure in under a minute. This script is what
# produced the confidence intervals, the paired tests, the calibration sweep and
# the cardinality result, and takes roughly **two hours on two cores**.
#
# Each script writes one row per run to results/, so the analysis in analysis/
# can be redone without re-running any simulation.
#
set -euo pipefail

cd "$(dirname "$0")"

PY="${PYTHON:-python3}"

echo "=============================================================="
echo " yamanaka-replay -- large-scale campaign (hours, not minutes)"
echo "=============================================================="
"$PY" --version
"$PY" -c "import numpy; print('numpy', numpy.__version__)"

mkdir -p results logs

START=$SECONDS

"$PY" analysis/test_stats.py            # the statistics layer checks itself first

"$PY" experiments/05_scale_main.py        # 1,200 runs -- intervals, paired tests
"$PY" experiments/06_synergy_threshold.py # 1,100 runs -- synergy vs calibration
"$PY" experiments/09_prior_quality.py     # 1,800 runs -- prior strength grid
"$PY" experiments/07_stopping_rule.py     # 1,680 runs -- stopping and calibration
"$PY" experiments/08_scaling.py           #   300 runs -- 24 to 60 candidates
"$PY" experiments/10_unknown_cardinality.py  # 240 runs -- the cardinality fix

"$PY" analysis/analyze.py | tee results/analysis_report.txt

echo "=============================================================="
echo " done in $((SECONDS - START))s -- per-run CSVs in results/"
echo "=============================================================="
