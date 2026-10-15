"""
The gap response as a latent variable, inferred alongside the recipe.

Why this module exists
----------------------
Experiment 13 measured what happens when the designer's assumed gap response
``rho`` is wrong in the way measured data says it is wrong: 2.5% of recipes
recovered, and the failures arrive with a mean posterior confidence of 0.991.
Told the true ``rho``, the same designer returns to 99%. So the published
method is sound and its calibration is not, which makes estimating ``rho`` a
precondition for using the method at all rather than a refinement of it.

There are two ways to meet that precondition.

**Calibrate first.** Measure ``rho`` on a system already solved -- a complete
recipe and each of its one-short pools -- then fix it and search. This is what
a bench scientist would do, and :func:`pilot_estimate_rho` implements it.

**Infer it in the loop.** Carry ``rho`` in the posterior with the recipe and
let the same experiments inform both. That is this module.

The second is only attractive if it is not much more expensive than the first,
and that is the open question experiment 16 answers. It is not obvious either
way: joint inference spends evidence identifying a nuisance parameter instead
of the recipe, while the pilot spends whole experiments up front and then
commits to a point estimate carrying its own sampling error.

The collapse still holds, one index wider
-----------------------------------------
The fixed designer is cheap because the colony distribution depends on a
hypothesis only through the overlap ``|S n h|``, so both the likelihood and the
expected information gain reduce to a mixture over five overlap classes however
many hypotheses are alive. Adding a discretised gap model indexes that mixture
by one more coordinate: the classes become ``(overlap, model)``, and the cost
of choosing an experiment grows by the number of grid points rather than by the
size of the posterior. At the default grid that is 45 classes instead of 5.

Discretisation is a real approximation, not a presentational one. The true
``rho`` in experiment 16 is drawn from a continuous measured distribution and
will generally fall between grid points, so the designer remains mildly
misspecified even when inference works perfectly. That residual is part of what
the experiment measures.
"""

from __future__ import annotations

import numpy as np

from . import inference as inf
from .inference import HYPOTHESES, poisson_pmf
from .oracle import ColonyOracle, RECIPE_SIZE
from .strategies import structured_pools

#: Grid over the gap response. Spans the value this study assumed (0.10), the
#: range measured in Norman et al. (median 0.77, IQR 0.62-1.03), and the
#: sub-additive region above 1 where one constituent outperforms the pair,
#: which 29% of measured pairs occupy. Spacing is geometric below 1 because the
#: likelihood is far more sensitive to a factor-of-two error at 0.05 than at
#: 0.8.
DEFAULT_RHOS = (0.02, 0.05, 0.10, 0.20, 0.35, 0.55, 0.80, 1.00, 1.25)

#: Fall-off shapes the designer may entertain. Pair measurements constrain the
#: one-short step and say nothing beyond it, so a designer that fixes the shape
#: is making an untested assumption of the same kind that ``rho = 0.10`` was.
DEFAULT_SHAPES = ("gap",)


class GapGrid:
    """
    A discretised set of assay models, each a ``(shape, rho)`` pair.

    Parameters
    ----------
    rhos:
        Grid points for the gap response.
    shapes:
        ``"gap"`` (efficiency vanishes beyond one missing factor) and/or
        ``"geometric"`` (efficiency is ``rho`` per missing factor). Passing
        both makes the shape latent too.
    """

    def __init__(self, rhos=DEFAULT_RHOS, shapes=DEFAULT_SHAPES):
        self.models = [(s, float(r)) for s in shapes for r in rhos]
        self.rho = np.array([r for _, r in self.models])
        self.shape = [s for s, _ in self.models]

    def __len__(self) -> int:
        return len(self.models)

    def efficiency(self, missing: int, model: int) -> float:
        """Assay efficiency for a pool ``missing`` factors, under one model."""
        if missing <= 0:
            return 1.0
        shape, rho = self.models[model]
        if shape == "geometric":
            return float(rho**missing)
        return float(rho) if missing == 1 else 0.0

    def rate_table(self, pool_size: int) -> np.ndarray:
        """``(RECIPE_SIZE + 1, n_models)`` of assumed rates for this pool size."""
        penalty = float(np.exp(-ColonyOracle.DECAY * max(0, pool_size - RECIPE_SIZE)))
        return np.array(
            [
                [
                    ColonyOracle.BASE * self.efficiency(RECIPE_SIZE - o, m) * penalty
                    for m in range(len(self.models))
                ]
                for o in range(RECIPE_SIZE + 1)
            ]
        )

    def pmf_table(self, pool_size: int) -> np.ndarray:
        """``(RECIPE_SIZE + 1, n_models, n_counts)`` count distributions."""
        rates = self.rate_table(pool_size)
        return np.array(
            [[poisson_pmf(rates[o, m]) for m in range(rates.shape[1])]
             for o in range(rates.shape[0])]
        )


# ------------------------------------------------------------ joint posterior


def joint_prior(grid: GapGrid, recipe_prior: np.ndarray | None = None) -> np.ndarray:
    """
    ``(n_hypotheses, n_models)`` prior, uniform over models by default.

    A uniform prior over the grid is deliberately uninformative about ``rho``:
    the point of the exercise is to see what the experiments alone can
    establish, and seeding the grid with the measured distribution would be
    assuming most of the answer.
    """
    n_h = len(HYPOTHESES)
    recipe = np.full(n_h, 1.0 / n_h) if recipe_prior is None else recipe_prior
    joint = np.outer(recipe, np.full(len(grid), 1.0 / len(grid)))
    return joint / joint.sum()


def _class_mass(posterior: np.ndarray, overlap: np.ndarray) -> np.ndarray:
    """``(RECIPE_SIZE + 1, n_models)`` posterior mass in each class."""
    return np.array(
        [posterior[overlap == o].sum(axis=0) for o in range(RECIPE_SIZE + 1)]
    )


def joint_update(
    grid: GapGrid, posterior: np.ndarray, pool, observed: int
) -> np.ndarray:
    """Posterior over ``(recipe, model)`` after seeing ``observed`` colonies."""
    overlap = inf._overlaps(pool)
    pmf = grid.pmf_table(len(pool))
    likelihood = pmf[:, :, min(observed, pmf.shape[2] - 1)]  # (overlaps, models)
    updated = posterior * likelihood[overlap]
    total = updated.sum()
    if total <= 0:
        # No model in the grid can explain the observation. Resetting to the
        # prior is the honest response: silently keeping a posterior that the
        # data has excluded is how the fixed designer produced 0.991
        # confidence on wrong answers.
        return joint_prior(grid)
    return updated / total


def joint_eig(grid: GapGrid, posterior: np.ndarray, pool) -> float:
    """
    Mutual information in bits between ``(recipe, model)`` and the outcome.

    Computed over the ``(overlap, model)`` classes, so the cost depends on the
    grid size and not on how concentrated the posterior is.
    """
    overlap = inf._overlaps(pool)
    mass = _class_mass(posterior, overlap)
    if mass.sum() <= 0:
        return 0.0
    pmf = grid.pmf_table(len(pool))

    flat_mass = mass.reshape(-1)
    flat_pmf = pmf.reshape(-1, pmf.shape[2])
    marginal = flat_mass @ flat_pmf
    conditional = float(
        sum(w * inf._entropy_bits(p) for w, p in zip(flat_mass, flat_pmf) if w > 0)
    )
    return inf._entropy_bits(marginal) - conditional


def marginal_recipe(posterior: np.ndarray) -> np.ndarray:
    """``P(recipe)``, with the assay model integrated out."""
    return posterior.sum(axis=1)


def marginal_model(posterior: np.ndarray) -> np.ndarray:
    """``P(model)``, with the recipe integrated out."""
    return posterior.sum(axis=0)


def posterior_rho(grid: GapGrid, posterior: np.ndarray) -> tuple[float, float, float]:
    """
    Posterior mean of ``rho`` and an equal-tailed 90% credible interval.

    The interval runs between *cell boundaries*, not between grid points. A
    discretised posterior places its mass on nine exact values while the
    quantity being estimated is continuous, so an interval reported as
    ``[0.80, 0.80]`` would exclude a true 0.77 by construction and the coverage
    of such an interval measures the spacing of the grid rather than the
    quality of the inference. Each grid point is therefore read as standing for
    the interval midway to its neighbours, which is what the discretisation
    means, with the outermost cells extended by half their adjacent spacing.

    Several models can share a ``rho`` when the shape is latent as well; their
    mass is pooled, since the shape is a separate question.
    """
    weights = marginal_model(posterior)
    values = np.unique(grid.rho)
    mass = np.array([weights[grid.rho == v].sum() for v in values])
    total = max(mass.sum(), 1e-12)
    mass = mass / total
    mean = float(values @ mass)

    midpoints = (values[:-1] + values[1:]) / 2.0
    lower_edges = np.concatenate(([values[0] - (midpoints[0] - values[0])], midpoints))
    upper_edges = np.concatenate((midpoints, [values[-1] + (values[-1] - midpoints[-1])]))

    cdf = np.cumsum(mass)
    first = int(np.searchsorted(cdf, 0.05))
    last = min(int(np.searchsorted(cdf, 0.95)), values.size - 1)
    return mean, float(lower_edges[first]), float(upper_edges[last])


# ------------------------------------------------------------------ designer


def adaptive_design_latent_gap(
    oracle,
    grid: GapGrid | None = None,
    seed: int = 0,
    max_experiments: int = 120,
    confidence: float = 0.95,
    recipe_prior: np.ndarray | None = None,
    n_pools: int = 120,
) -> tuple[set[int], int, float, np.ndarray, GapGrid]:
    """
    Information-gain design over a joint posterior on recipe and assay model.

    Stopping is on the **recipe** marginal, not the joint: ``rho`` is a
    nuisance parameter, and a designer that waited to pin it down would spend
    experiments answering a question nobody asked.

    Returns ``(recipe, experiments, recipe confidence, posterior, grid)``.
    """
    grid = grid or GapGrid()
    rng = np.random.default_rng(seed)
    posterior = joint_prior(grid, recipe_prior)

    while oracle.n_experiments < max_experiments:
        recipe = marginal_recipe(posterior)
        if recipe.max() >= confidence:
            break
        pools = structured_pools(rng, recipe, n_pools)
        pool = max(pools, key=lambda p: joint_eig(grid, posterior, p))
        posterior = joint_update(grid, posterior, pool, oracle.run(pool))

    recipe = marginal_recipe(posterior)
    best = HYPOTHESES[int(np.argmax(recipe))]
    return set(best), oracle.n_experiments, float(recipe.max()), posterior, grid


# ------------------------------------------------------------------- pilot


def pilot_estimate_rho(oracle, reference, replicates: int = 1) -> tuple[float, int]:
    """
    Calibrate ``rho`` the way a laboratory would, on a recipe already known.

    Runs the complete reference recipe and each of its one-short pools, and
    returns the ratio of mean one-short yield to full yield, with the
    experiments it consumed.

    **This arm is handed something the latent designer is not.** In this
    benchmark the only recipe that yields colonies is the hidden one, so the
    pilot must be pointed at it; a real calibration would instead use a
    different, previously solved system. The arm is therefore an *upper bound*
    on what a pilot can deliver -- it pays the experiments but not the risk of
    calibrating on the wrong system. If joint inference matches it anyway, the
    comparison is settled in the harder direction.

    Returns ``(rho estimate, experiments used)``.
    """
    reference = sorted(reference)
    start = oracle.n_experiments
    full = [oracle.run(reference) for _ in range(replicates)]
    short = [
        oracle.run([g for g in reference if g != drop])
        for drop in reference
        for _ in range(replicates)
    ]
    denominator = float(np.mean(full))
    if denominator <= 0:
        # The positive control gave nothing: the pilot has learned only that,
        # and returning the assumed value is what a lab with no signal would
        # fall back on.
        return float(ColonyOracle.PARTIAL_3), oracle.n_experiments - start
    return float(np.mean(short)) / denominator, oracle.n_experiments - start
