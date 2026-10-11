"""
Colony-formation oracle: a calibrated stand-in for Yamanaka's selection assay.

Takahashi & Yamanaka (Cell, 2006) did not optimise a continuous score. Their
readout was the number of drug-resistant colonies surviving Fbx15/G418
selection: an *integer*, *stochastic*, and *very expensive* signal, where a
single evaluation is a wet-lab experiment taking weeks.

This module reproduces that readout as a Poisson process whose rate is set by
three effects documented in the original paper:

1.  **All-or-nothing core.** A pool containing all four essential factors
    reprograms efficiently. A pool missing any of them does not.
2.  **Graded partial effect.** Pools containing exactly three of the four gave
    a small number of colonies that lacked full ES-like morphology. This
    partial signal is what makes the search space navigable at all; the
    `03_misspecification.py` experiment shows the method collapses without it.
3.  **Pool-size penalty.** Transducing more factors lowers per-cell efficiency.
    This is why the paper's 10-factor subset outperformed the full 24.

Calibration targets taken from the paper:

    24 factors -> ~22 colonies   (figure reported in the paper)
    10 factors -> ~36 colonies   ("the 10 outperformed the 24")
     4 factors -> ~45 colonies   (clean recipe)
     1 factor  ->   0 colonies   (reported: no single factor sufficed)
     3 of 4    ->  ~4 colonies   (partial effect)
"""

from __future__ import annotations

import numpy as np

N_CANDIDATES = 24  # Yamanaka's shortlist
RECIPE_SIZE = 4  # OCT4 + SOX2 + KLF4 + MYC


class ColonyOracle:
    """
    The laboratory. Returns a colony count for a transduced factor pool.

    The expected count is

        lambda(S) = BASE * f(|S n R|) * exp(-DECAY * max(0, |S| - 4))

    where ``R`` is the hidden essential set, ``f(4) = 1``, ``f(3) = PARTIAL_3``
    and ``f(k) = 0`` for ``k < 3``. Observed counts are Poisson draws around
    that rate, so repeated identical experiments disagree, exactly as bench
    replicates do.

    Parameters
    ----------
    truth:
        The hidden essential factor set. Searchers never see this.
    seed:
        Seed for the observation noise.

    Attributes
    ----------
    n_experiments:
        Count of wet-lab experiments consumed. **This is the metric the whole
        study optimises** -- not wall-clock time, not FLOPs.
    log:
        ``(pool, observed_count)`` for every experiment run, in order.
    """

    BASE = 45.0
    DECAY = 0.0357  # calibrated so that a 24-factor pool yields ~22 colonies
    PARTIAL_3 = 0.10

    def __init__(self, truth, seed: int = 0):
        self.truth = frozenset(truth)
        self._rng = np.random.default_rng(seed)
        self.n_experiments = 0
        self.log: list[tuple[tuple[int, ...], int]] = []

    def rate(self, pool) -> float:
        """Expected colony count for ``pool`` (no noise). Used for calibration."""
        pool = set(pool)
        if not pool:
            return 0.0
        hits = len(pool & self.truth)
        if hits == RECIPE_SIZE:
            efficiency = 1.0
        elif hits == RECIPE_SIZE - 1:
            efficiency = self.PARTIAL_3
        else:
            efficiency = 0.0
        size_penalty = np.exp(-self.DECAY * max(0, len(pool) - RECIPE_SIZE))
        return self.BASE * efficiency * size_penalty

    def run(self, pool) -> int:
        """Run one wet-lab experiment on ``pool`` and return the colony count."""
        self.n_experiments += 1
        count = int(self._rng.poisson(self.rate(pool)))
        self.log.append((tuple(sorted(pool)), count))
        return count


class MisspecifiedOracle(ColonyOracle):
    """
    An oracle whose biology deliberately disagrees with the searcher's model.

    Used by ``03_misspecification.py`` to find the boundary of the method.
    Every parameter left as ``None`` keeps the calibrated default.

    Parameters
    ----------
    partial_3, base, decay:
        Override the corresponding oracle constants.
    extra_required:
        A fifth factor that must also be present for any colonies to form --
        i.e. the true recipe is larger than the searcher assumes. This is the
        one perturbation the method cannot survive.
    """

    def __init__(
        self,
        truth,
        seed: int = 0,
        partial_3: float | None = None,
        base: float | None = None,
        decay: float | None = None,
        extra_required: int | None = None,
    ):
        super().__init__(truth, seed)
        if partial_3 is not None:
            self.PARTIAL_3 = partial_3
        if base is not None:
            self.BASE = base
        if decay is not None:
            self.DECAY = decay
        self.extra_required = extra_required

    def rate(self, pool) -> float:
        if self.extra_required is not None and self.extra_required not in set(pool):
            return 0.0
        return super().rate(pool)
