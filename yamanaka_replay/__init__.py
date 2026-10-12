"""
yamanaka-replay
===============

Replaying the 2006 Yamanaka-factor discovery in silico, under the original
colony-formation fitness function, to measure how many *wet-lab experiments*
adaptive experimental design would have needed.

Headline: a median of 6-11 experiments against 61 for the historical
leave-one-out protocol.

The oracle is a simulator calibrated to published figures, not an experiment.
No biological claim is made or supported by this code.
"""

from .oracle import ColonyOracle, MisspecifiedOracle, N_CANDIDATES, RECIPE_SIZE
from .inference import (
    HYPOTHESES,
    bayes_update,
    expected_information_gain,
    marginal_factor_probability,
    uniform_prior,
    weighted_prior,
)
from .strategies import adaptive_design, historical_protocol

__version__ = "1.1.0"

__all__ = [
    "ColonyOracle",
    "MisspecifiedOracle",
    "N_CANDIDATES",
    "RECIPE_SIZE",
    "HYPOTHESES",
    "bayes_update",
    "expected_information_gain",
    "marginal_factor_probability",
    "uniform_prior",
    "weighted_prior",
    "adaptive_design",
    "historical_protocol",
]
