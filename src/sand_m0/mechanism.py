"""Exact deterministic structure of M0 (derived in docs/M0_MECHANISM.md).

For the M0 update rule with a fixed near zone (x < a) and respray shape q,
the expected-mass dynamics collapse onto a one-dimensional path. Writing
Q = (near-zone mass of q), Q_f = 1 - Q, and

    z_{t+1} = rho * z_t,   rho = 1 - alpha * (1 - Q),   z_0 = 1,

the exact solution (proved by induction) is

    near bins:  m_i(t) = z_t * q_i
    far bins:   m_i(t) = lambda_t * q_i,
                lambda_t = (1 - Q * z_t) / Q_f   (M0-A: 4 - 3 z_t)

so the whole deterministic trajectory is a function of the single state
variable z ("how much of the original near-zone mass is left"). The
histogram uniformity error along the path is an exact quadratic,

    U(z) = S_n * z^2 + S_f * ((1 - Q z) / Q_f)^2 - 1/K,

with S_n, S_f the exact bin-probability sums of squares over the two zones.
All functions here are exact algebra on the discrete bins — no numerical
approximation is introduced anywhere in this module.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import ThrowDistribution, near_zone_mask


@dataclass(frozen=True)
class MechanismConstants:
    """Exact zone constants of an M0 configuration (dist, K, a)."""

    K: int
    Q_near: float          # near-zone mass of q at t = 0
    Q_far: float           # far-zone mass of q
    S_near: float          # sum of q_i^2 over near bins (shape concentration)
    S_far: float           # sum of q_i^2 over far bins
    near_mask: np.ndarray
    q: np.ndarray          # bin probabilities of the throw distribution


def mechanism_constants(dist: ThrowDistribution, K: int, a: float) -> MechanismConstants:
    q = dist.bin_probs(K)
    near = near_zone_mask(K, a)
    return MechanismConstants(
        K=K,
        Q_near=float(q[near].sum()),
        Q_far=float(q[~near].sum()),
        S_near=float((q[near] ** 2).sum()),
        S_far=float((q[~near] ** 2).sum()),
        near_mask=near,
        q=q,
    )


def decay_factor(alpha: float, Q_near: float) -> float:
    """Per-round near-mass retention factor rho = 1 - alpha * (1 - Q_near).

    Each round the near zone loses an alpha fraction of its mass, but a
    Q_near share of every respray falls back into the near zone, so the net
    loss per round is alpha * (1 - Q_near).
    """
    return 1.0 - alpha * (1.0 - Q_near)


def near_amplitude(t: int, rho: float) -> float:
    """z_t = rho^t: scale of the near-zone shape at round t."""
    return rho ** t


def far_amplitude(z: float, Q_near: float, Q_far: float) -> float:
    """lambda(z) = (1 - Q_near * z) / Q_far — scale of the far-zone shape.

    Mass conservation: the far zone holds 1 - Q_near * z in total, spread
    over the fixed far-zone shape q, so the amplitude is that mass divided
    by Q_far. For M0-A (Q = 3/4, Q_f = 1/4) this is 4 - 3z.
    """
    return (1.0 - Q_near * z) / Q_far


def analytic_state(
    z: float, constants: MechanismConstants
) -> np.ndarray:
    """Exact expected bin-mass vector at the state labelled by z."""
    lam = far_amplitude(z, constants.Q_near, constants.Q_far)
    return np.where(constants.near_mask, z * constants.q, lam * constants.q)


def analytic_trajectory(
    dist: ThrowDistribution, *, K: int, T: int, a: float, alpha: float
) -> list[np.ndarray]:
    """Exact expected states for rounds t = 0..T (closed form, no stepping)."""
    c = mechanism_constants(dist, K, a)
    rho = decay_factor(alpha, c.Q_near)
    return [analytic_state(near_amplitude(t, rho), c) for t in range(T + 1)]


def quadratic_coefficients(constants: MechanismConstants) -> tuple[float, float, float]:
    """Coefficients (A, B, C) of the exact quadratic U(z) = A z^2 + B z + C.

    Expanded from U(z) = S_n z^2 + S_f ((1 - Q z)/Q_f)^2 - 1/K.
    """
    Qn, Qf, Sn, Sf, K = (
        constants.Q_near, constants.Q_far,
        constants.S_near, constants.S_far, constants.K,
    )
    A = Sn + (Qn / Qf) ** 2 * Sf
    B = -2.0 * Qn * Sf / Qf**2
    C = Sf / Qf**2 - 1.0 / K
    return A, B, C


def u_of_z(z: float, constants: MechanismConstants) -> float:
    A, B, C = quadratic_coefficients(constants)
    return A * z * z + B * z + C


def optimal_z(constants: MechanismConstants) -> float:
    """Minimizer of the exact quadratic: z* = Q S_f / (Q^2 S_f + Q_f^2 S_n)."""
    Qn, Qf, Sn, Sf = (
        constants.Q_near, constants.Q_far, constants.S_near, constants.S_far
    )
    return Qn * Sf / (Qn**2 * Sf + Qf**2 * Sn)


def optimal_round(rho: float, z_star: float) -> float:
    """Continuous best round t* = ln(z*) / ln(rho) (needs 0 < rho < 1)."""
    if not 0.0 < rho < 1.0:
        raise ValueError(f"t* undefined for rho={rho}: no time progression")
    import math

    return math.log(z_star) / math.log(rho)
