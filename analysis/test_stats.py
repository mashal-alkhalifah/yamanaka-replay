#!/usr/bin/env python3
"""
Checks on the statistics layer.

Every inferential claim in the analysis rests on these four functions, and they
are hand-rolled to keep the dependency list at one package. So each is checked
against a value that is known independently: a textbook figure, a closed-form
answer, or a case where the right answer is obvious by construction.

Run with ``python analysis/test_stats.py``; exits nonzero on failure.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from stats import (  # noqa: E402
    bootstrap_ci,
    logistic_threshold,
    loglog_slope,
    mcnemar_exact,
    paired_permutation_test,
    wilson_ci,
)

failures: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    print(f"  {'ok  ' if condition else 'FAIL'}  {name}{'   ' + detail if detail else ''}")
    if not condition:
        failures.append(name)


# ---- Wilson interval -------------------------------------------------------
# Checked against its definition rather than a quoted table, because quoted
# tables differ over whether a continuity correction is applied. The Wilson
# bounds are by construction the two roots p of
#       |p_hat - p| / sqrt(p (1 - p) / n) = z
# so each returned endpoint must drive that score statistic back to exactly z.
lo, hi = wilson_ci(15, 20)


def score_at(p: float, p_hat: float = 0.75, n: int = 20) -> float:
    return abs(p_hat - p) / np.sqrt(p * (1 - p) / n)


check("wilson_ci endpoints solve the defining score equation",
      abs(score_at(lo) - 1.959964) < 1e-3 and abs(score_at(hi) - 1.959964) < 1e-3,
      f"[{lo:.4f}, {hi:.4f}] -> scores {score_at(lo):.4f}, {score_at(hi):.4f}")

check("wilson_ci brackets the point estimate", lo < 0.75 < hi)

# At 100% success the interval must still have width -- the reason for using it.
lo, hi = wilson_ci(300, 300)
check("wilson_ci has nonzero width at 100% success", lo < 1.0 and hi == 1.0,
      f"[{lo:.3f}, {hi:.3f}]")

# ---- exact McNemar --------------------------------------------------------
# All discordant pairs on one side of a 10-pair split: 2 * 0.5^10 = 0.00195.
p = mcnemar_exact(both=50, only_a=10, only_b=0, neither=40)
check("mcnemar_exact reproduces 2*0.5^10 for a 10-0 discordant split",
      abs(p - 2 * 0.5**10) < 1e-9, f"p = {p:.6f}")

check("mcnemar_exact returns 1.0 with no discordant pairs",
      mcnemar_exact(10, 0, 0, 10) == 1.0)

check("mcnemar_exact is symmetric in its discordant counts",
      abs(mcnemar_exact(5, 3, 9, 5) - mcnemar_exact(5, 9, 3, 5)) < 1e-12)

# ---- paired permutation test ----------------------------------------------
rng = np.random.default_rng(0)
same = rng.normal(size=400)
_, p_null = paired_permutation_test(same, same + rng.normal(scale=1e-9, size=400))
check("permutation test does not reject when the samples are identical",
      p_null > 0.05, f"p = {p_null:.3f}")

shifted = same + 2.0
diff, p_alt = paired_permutation_test(shifted, same)
check("permutation test recovers a known +2.0 shift and rejects",
      abs(diff - 2.0) < 0.01 and p_alt < 1e-3, f"diff = {diff:.3f}, p = {p_alt:.5f}")

# False-positive rate under the null should sit near the nominal 5%.
rejections = 0
for trial in range(200):
    r = np.random.default_rng(1000 + trial)
    x, y = r.normal(size=60), r.normal(size=60)
    if paired_permutation_test(x, y, n_resamples=2000, seed=trial)[1] < 0.05:
        rejections += 1
check("permutation test holds its nominal size under the null",
      0.01 <= rejections / 200 <= 0.11, f"{rejections / 200:.3f} rejected at alpha=0.05")

# ---- bootstrap ------------------------------------------------------------
# Median of a large standard normal sample: the interval must cover 0.
sample = np.random.default_rng(3).normal(size=2000)
lo, hi = bootstrap_ci(sample, np.median, n_resamples=4000)
check("bootstrap median interval covers the true median of 0",
      lo < 0 < hi, f"[{lo:.3f}, {hi:.3f}]")

# Width must shrink roughly as 1/sqrt(n).
w_small = np.subtract(*reversed(bootstrap_ci(sample[:200], np.mean, n_resamples=4000)))
w_large = np.subtract(*reversed(bootstrap_ci(sample, np.mean, n_resamples=4000)))
check("bootstrap interval narrows with sample size",
      w_large < w_small, f"{w_small:.3f} -> {w_large:.3f}")

# ---- logistic threshold ---------------------------------------------------
# Generate data from a known logistic in log10(x) and recover its 50% point.
true_threshold = 0.03
x_grid = np.repeat(np.array([0.003, 0.01, 0.02, 0.03, 0.05, 0.1, 0.3]), 400)
slope_true = 4.0
eta = slope_true * (np.log10(x_grid) - np.log10(true_threshold))
y = (np.random.default_rng(5).random(x_grid.size) < 1 / (1 + np.exp(-eta))).astype(float)
est, slope = logistic_threshold(x_grid, y, 0.5)
check("logistic_threshold recovers a planted 50% crossing point",
      abs(est - true_threshold) / true_threshold < 0.15, f"{est:.4f} vs {true_threshold}")
check("logistic_threshold recovers the planted slope sign and scale",
      abs(slope - slope_true) / slope_true < 0.25, f"{slope:.2f} vs {slope_true}")

check("logistic_threshold refuses degenerate input",
      np.isnan(logistic_threshold(np.array([1.0, 2.0]), np.array([1.0, 1.0]))[0]))

# ---- log-log slope --------------------------------------------------------
x = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
check("loglog_slope recovers a quadratic exponent", abs(loglog_slope(x, x**2) - 2.0) < 1e-9)
check("loglog_slope recovers a square-root exponent",
      abs(loglog_slope(x, np.sqrt(x)) - 0.5) < 1e-9)

print()
if failures:
    print(f"  {len(failures)} check(s) failed: {', '.join(failures)}")
    sys.exit(1)
print("  all statistics checks passed")
