"""
Oracles grounded in measured data rather than in the 2006 paper's constants.

Why this module exists
----------------------
The simulated :class:`~yamanaka_replay.oracle.ColonyOracle` is calibrated from
Takahashi & Yamanaka (2006) itself, so replaying that paper against it cannot
validate it -- the target was used to set the parameters. The audit in
``analysis/real_data_audit.py`` then showed that the one assumption the whole
study rests on is the one measured data refuses: the oracle assumes a pool one
factor short of the recipe performs at ``rho = 0.10`` of a complete pool, while
in Norman et al. (2019) the measured analogue has median 0.77 and never falls
below 0.20 in 131 combinations.

This module closes that loop in the two ways the data allows.

``EmpiricalGapOracle``
    Keeps the paper's task -- find a hidden recipe among ``n`` candidates using
    as few assays as possible -- but takes the gap response ``rho`` from the
    measured Norman distribution instead of from the 2006 calibration. Each run
    draws a different ``rho``, so a campaign of runs reports performance over
    the landscape biology actually presents rather than at one assumed point.

``NormanOracle``
    Spends no simulation at all. Every answer it returns is a measured
    transcriptome from Norman et al. (2019). It will only answer queries for
    conditions that dataset actually contains, which restricts pools to the
    singles and pairs that were profiled, and so restricts the task to
    discovery within a measured library rather than open-ended pool design.
    That restriction is the price of using only real numbers, and it is
    stated rather than engineered around.

The identification assumption
-----------------------------
Norman measures transcriptional response magnitude; the oracle speaks in colony
counts. Carrying ``rho`` from one to the other assumes response magnitude is
monotone in assay yield. That assumption is not tested here and cannot be
tested with this data -- it is the hinge on which every number in
``EmpiricalGapOracle`` turns, and it is stated here so a reader can discount
the results accordingly. ``NormanOracle`` does not need it, which is the main
reason it is worth the loss of generality.
"""

from __future__ import annotations

import csv
import pathlib

import numpy as np

from .oracle import ColonyOracle, RECIPE_SIZE

_REPO = pathlib.Path(__file__).resolve().parents[1]
#: One-short ratios measured on Norman pairs, written by ``real_data_audit.py``.
AUDIT_CSV = _REPO / "results" / "real_data_audit_norman.csv"

#: Relative discrepancy between the two halves of the gene panel, used as the
#: measurement noise scale for :class:`NormanOracle`. Recomputed from the data
#: by :func:`estimate_split_half_cv`, which ``experiments/14`` prints so the
#: constant can be checked against the dataset it came from.
SPLIT_HALF_CV = 0.05


def estimate_split_half_cv(pseudobulk: pathlib.Path, seed: int = 0) -> float:
    """
    Measurement noise, taken from the dataset's disagreement with itself.

    The gene panel is split at random into two halves and each condition is
    scored on each half alone, rescaled to be comparable with the full-panel
    score. The median relative discrepancy between the two is returned. This is
    a lower bound on true assay noise -- two halves of one library share
    everything upstream of the measurement -- but it has the virtue of being
    measured rather than chosen, and a noise level chosen by the person
    reporting the result is worth very little.
    """
    labels, matrix = [], []
    with pathlib.Path(pseudobulk).open(newline="") as fh:
        reader = csv.reader(fh)
        next(reader)
        for row in reader:
            labels.append(row[0])
            matrix.append([float(v) for v in row[1:]])
    X = np.asarray(matrix, dtype=float)
    control = next(
        (i for i, lab in enumerate(labels) if lab in ("control", "ctrl", "non-targeting")),
        None,
    )
    if control is None:
        raise ValueError("no control condition in the pseudobulk table")
    delta = X - X[control]

    order = np.random.default_rng(seed).permutation(delta.shape[1])
    half = delta.shape[1] // 2
    scale = np.sqrt(delta.shape[1] / half)
    a = np.linalg.norm(delta[:, order[:half]], axis=1) * scale
    b = np.linalg.norm(delta[:, order[half:]], axis=1) * scale
    full = np.maximum(np.linalg.norm(delta, axis=1), 1e-9)
    return float(np.median(np.abs(a - b) / full))


# ------------------------------------------------------- measured gap profile


def measured_one_short(path: pathlib.Path | None = None) -> np.ndarray:
    """
    The measured one-short ratios, as a plain array.

    These are ``max_constituent ||delta|| / ||delta_combination||`` over the 131
    Norman pairs whose singles were also profiled. Values at or above 1 are kept
    rather than clipped: they mean one constituent moved the cell as far as the
    pair did, and a design method that assumes otherwise should be made to face
    them.
    """
    path = path or AUDIT_CSV
    with path.open(newline="") as fh:
        values = [
            float(row["one_short"])
            for row in csv.DictReader(fh)
            if row.get("one_short") not in (None, "")
        ]
    if not values:
        raise ValueError(f"no one-short values in {path}")
    return np.asarray(values, dtype=float)


class EmpiricalGapOracle(ColonyOracle):
    """
    The 2006 assay with its gap response replaced by a measured one.

    Parameters
    ----------
    truth:
        The hidden recipe.
    seed:
        Seed for both the draw of ``rho`` and the Poisson observation noise.
    rho:
        Fix the gap response instead of drawing it. Used to hold ``rho``
        constant while something else varies.
    decay:
        How efficiency falls as a pool loses further factors.

        ``"gap"``
            The paper's shape: a pool three-of-four gives ``rho``, anything
            less gives nothing. Preserves the sharp threshold and changes only
            its depth.
        ``"geometric"``
            ``efficiency = rho ** (recipe_size - overlap)``, so a pool two
            short still gives ``rho**2``. At the measured ``rho`` this is a far
            flatter landscape, with signal everywhere and therefore much weaker
            evidence in any single assay. Which shape is right is not something
            pair data can settle, so both are reported.

    Attributes
    ----------
    rho:
        The gap response this instance drew, recorded so analyses can condition
        on it.
    """

    def __init__(
        self,
        truth,
        seed: int = 0,
        rho: float | None = None,
        decay: str = "gap",
        one_short: np.ndarray | None = None,
    ):
        super().__init__(truth, seed)
        if decay not in ("gap", "geometric"):
            raise ValueError("decay must be 'gap' or 'geometric'")
        self.decay_shape = decay
        if rho is None:
            sample = measured_one_short() if one_short is None else one_short
            rho = float(np.random.default_rng(seed).choice(sample))
        self.rho = float(rho)

    def rate(self, pool) -> float:
        pool = set(pool)
        if not pool:
            return 0.0
        hits = len(pool & self.truth)
        short_by = RECIPE_SIZE - hits
        if short_by == 0:
            efficiency = 1.0
        elif self.decay_shape == "gap":
            efficiency = self.rho if short_by == 1 else 0.0
        else:
            efficiency = self.rho**short_by
        penalty = np.exp(-self.DECAY * max(0, len(pool) - RECIPE_SIZE))
        return self.BASE * efficiency * penalty


# ------------------------------------------------------- strictly measured


class NormanOracle:
    """
    A retrospective assay: every answer is a measurement, none is simulated.

    The phenotype returned for a condition is the magnitude of its pseudobulk
    expression shift against the unperturbed control, ``||x - x_control||``.
    That is a deliberately plain choice. Any score tuned to reward a particular
    biology would make the search look good by construction, and the point of a
    retrospective oracle is to remove exactly that freedom.

    Only conditions present in the dataset can be queried. Asking for anything
    else raises :class:`KeyError` rather than interpolating, because an
    interpolated answer would reintroduce the simulation this class exists to
    avoid.

    Observation noise is multiplicative lognormal at :data:`SPLIT_HALF_CV`,
    estimated by scoring each condition on two disjoint halves of the gene panel
    and comparing. It is the dataset's own internal disagreement, not a figure
    chosen to make the problem interesting.

    Parameters
    ----------
    pseudobulk, perturbations:
        Paths to the files written by ``data/prepare_norman.py``.
    seed:
        Seed for the observation noise.
    noiseless:
        Return the measured value directly. For computing ground truth, never
        for a search -- a searcher given noiseless access is not solving the
        problem the wet lab poses.
    """

    def __init__(
        self,
        pseudobulk: pathlib.Path,
        perturbations: pathlib.Path,
        seed: int = 0,
        noiseless: bool = False,
    ):
        labels, matrix = [], []
        with pathlib.Path(pseudobulk).open(newline="") as fh:
            reader = csv.reader(fh)
            next(reader)
            for row in reader:
                labels.append(row[0])
                matrix.append([float(v) for v in row[1:]])
        X = np.asarray(matrix, dtype=float)
        index = {label: i for i, label in enumerate(labels)}

        control = next(
            (c for c in ("control", "ctrl", "non-targeting") if c in index), None
        )
        if control is None:
            raise ValueError("no control condition in the pseudobulk table")

        rows = list(csv.DictReader(pathlib.Path(perturbations).open(newline="")))
        self.genes_of: dict[frozenset[str], str] = {}
        for row in rows:
            if row["type"] == "control" or row["perturbation"] not in index:
                continue
            self.genes_of[frozenset(row["genes"].split(";"))] = row["perturbation"]

        delta = X - X[index[control]]
        self._score = {
            key: float(np.linalg.norm(delta[index[label]]))
            for key, label in self.genes_of.items()
        }

        self.candidates: list[str] = sorted(
            {g for key in self.genes_of for g in key}
        )
        self._rng = np.random.default_rng(seed)
        self.noiseless = noiseless
        self.n_experiments = 0
        self.log: list[tuple[tuple[str, ...], float]] = []

    # -- queries

    def available(self) -> list[frozenset[str]]:
        """Every condition that can be queried, as gene sets."""
        return list(self._score)

    def truth(self, pool) -> float:
        """The measured value with no observation noise. Not for searchers."""
        return self._score[frozenset(pool)]

    def best(self) -> tuple[frozenset[str], float]:
        """The strongest measured condition and its value -- the target."""
        key = max(self._score, key=self._score.get)
        return key, self._score[key]

    def run(self, pool) -> float:
        """
        Spend one experiment on ``pool`` and return its measured phenotype.

        Raises
        ------
        KeyError
            If that condition was never profiled. The caller must choose among
            :meth:`available` conditions; this oracle does not invent data.
        """
        key = frozenset(pool)
        if key not in self._score:
            raise KeyError(f"{sorted(key)} was not measured in this dataset")
        self.n_experiments += 1
        value = self._score[key]
        if not self.noiseless:
            value *= float(np.exp(self._rng.normal(0.0, SPLIT_HALF_CV)))
        self.log.append((tuple(sorted(key)), value))
        return value
