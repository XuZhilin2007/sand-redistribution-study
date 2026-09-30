"""Canonical M0 model: one-dimensional selective recycle-and-respray process.

Model version: M0 (first canonical baseline of the sand redistribution study,
established 2026-09-17). See docs/M0_MODEL.md for the full Chinese model
definition and evidence classification.

Space
-----
x in [0, 1] is an equal-area coordinate: equal-length intervals of x cover
equal ground area, so a perfectly spread state corresponds to the uniform
distribution on [0, 1]. x = 0 is nearest to the operator, x = 1 is farthest.

Update rule (fixed, applied once per round)
-------------------------------------------
1. Let M_near(t) be the total mass with x < a (the near zone).
2. Each unit of near-zone mass is selected independently with probability
   alpha and removed.
3. Every removed unit is resampled independently from the experiment's throw
   distribution q(x) on [0, 1].
4. Unselected mass stays where it is.

The M0 baseline parameters (a = 0.50, alpha = 0.25, K = 100, T = 30) are
toy-model choices fixed in experiments/m0_baseline/config.json. They are not
real-world optima and must not be described as such.
"""

__version__ = "0.1.0"

from .adaptive import (
    cumulative_excess,
    deficit_fill_redistribution,
    run_m1a_deterministic,
    run_m1a_mc,
)
from .kernel_theory import (
    a_crit_gate,
    contraction_margin,
    cumulative_excess_,
    cumulative_mismatch,
    mismatch_reinjection,
    one_step_excess,
    prefix_cdf,
    rebound_criterion,
    rebound_residual,
    state_reduced_gate,
    target_matched_excess,
)
from .diagnostics import (
    boundary_diagnostics,
    largest_jumps,
    recurrence_l1,
    sup_cdf_minus_uniform,
    window_stats,
)
from .mechanism import (
    MechanismConstants,
    analytic_state,
    analytic_trajectory,
    decay_factor,
    far_amplitude,
    mechanism_constants,
    near_amplitude,
    optimal_round,
    optimal_z,
    quadratic_coefficients,
    u_of_z,
)
from .model import (
    DISTRIBUTIONS,
    ThrowDistribution,
    near_zone_mask,
    noise_floor_l2,
    power_distribution,
    total_variation,
    uniformity_l2,
)
from .simulate import run_deterministic, run_deterministic_states, run_mc

__all__ = [
    "DISTRIBUTIONS",
    "ThrowDistribution",
    "cumulative_excess",
    "deficit_fill_redistribution",
    "run_m1a_deterministic",
    "run_m1a_mc",
    "mismatch_reinjection",
    "one_step_excess",
    "prefix_cdf",
    "target_matched_excess",
    "a_crit_gate",
    "contraction_margin",
    "cumulative_excess_",
    "cumulative_mismatch",
    "rebound_criterion",
    "rebound_residual",
    "state_reduced_gate",
    "boundary_diagnostics",
    "largest_jumps",
    "recurrence_l1",
    "sup_cdf_minus_uniform",
    "window_stats",
    "MechanismConstants",
    "analytic_state",
    "analytic_trajectory",
    "decay_factor",
    "far_amplitude",
    "mechanism_constants",
    "near_amplitude",
    "optimal_round",
    "optimal_z",
    "quadratic_coefficients",
    "u_of_z",
    "near_zone_mask",
    "noise_floor_l2",
    "power_distribution",
    "total_variation",
    "uniformity_l2",
    "run_deterministic",
    "run_deterministic_states",
    "run_mc",
    "__version__",
]
