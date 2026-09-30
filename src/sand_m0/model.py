"""Core M0 definitions: throw distributions, near zone, uniformity metrics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class ThrowDistribution:
    """A throw distribution q(x) on [0, 1] with exact CDF and inverse CDF.

    The inverse CDF is used for exact sampling (no lookup tables), and the
    CDF is differenced on bin edges to get exact bin probabilities for the
    deterministic expected-mass version.
    """

    name: str
    label: str
    pdf: Callable[[np.ndarray], np.ndarray]
    cdf: Callable[[np.ndarray], np.ndarray]
    ppf: Callable[[np.ndarray], np.ndarray]

    def bin_probs(self, K: int) -> np.ndarray:
        """Exact probability mass of each of K equal-width bins on [0, 1]."""
        edges = np.linspace(0.0, 1.0, K + 1)
        probs = np.diff(self.cdf(edges))
        # CDF differences telescope to 1; correct floating-point residue only.
        probs = probs / probs.sum()
        return probs

    def sample(self, rng: np.random.Generator, n: int) -> np.ndarray:
        """Draw n independent samples from q(x) via the inverse CDF."""
        return self.ppf(rng.random(n))


def _near_cdf(x: np.ndarray) -> np.ndarray:
    return 2.0 * x - x * x


def _near_ppf(u: np.ndarray) -> np.ndarray:
    return 1.0 - np.sqrt(1.0 - u)


def _far_cdf(x: np.ndarray) -> np.ndarray:
    return x * x


def _far_ppf(u: np.ndarray) -> np.ndarray:
    return np.sqrt(u)


DISTRIBUTIONS = {
    "near": ThrowDistribution(
        name="near",
        label="q_near(x) = 2(1-x) (near-biased)",
        pdf=lambda x: 2.0 * (1.0 - x),
        cdf=_near_cdf,
        ppf=_near_ppf,
    ),
    "uniform": ThrowDistribution(
        name="uniform",
        label="q_uniform(x) = 1 (uniform control)",
        pdf=lambda x: np.ones_like(x),
        cdf=lambda x: x,
        ppf=lambda u: u,
    ),
    "far": ThrowDistribution(
        name="far",
        label="q_far(x) = 2x (far-biased control)",
        pdf=lambda x: 2.0 * x,
        cdf=_far_cdf,
        ppf=_far_ppf,
    ),
}


def power_distribution(p: float) -> ThrowDistribution:
    """Member of the Phase 2A single-parameter power family.

        q_p(x) = (p + 1) * (1 - x)^p,   p > -1

    p = 0 is the uniform distribution, p = 1 is exactly the M0 baseline
    q_near, p > 0 tilts toward the operator (near-biased), -1 < p < 0
    tilts away (far-biased). The parameter is the power-law tilt of the
    density toward x = 0. Closed forms require p > -1 for normalisation
    and p > -0.5 for a finite continuous integral of q^2 (the discrete
    implementation only needs p > -1).

    Exact CDF (F(x) = 1 - (1-x)^(p+1)) and inverse CDF
    (x = 1 - (1-u)^(1/(p+1))) are analytic, so sampling and bin
    probabilities are stable for the whole family.
    """
    p = float(p)
    if p <= -1.0:
        raise ValueError(f"power family requires p > -1, got {p}")
    one_plus_p = p + 1.0

    def cdf(x: np.ndarray) -> np.ndarray:
        return 1.0 - (1.0 - np.asarray(x, dtype=float)) ** one_plus_p

    def ppf(u: np.ndarray) -> np.ndarray:
        return 1.0 - (1.0 - np.asarray(u, dtype=float)) ** (1.0 / one_plus_p)

    return ThrowDistribution(
        name=f"power_p{p:g}",
        label=f"q_p(x) = {p + 1:g}(1-x)^{p:g} (power family, p={p:g})",
        pdf=None,
        cdf=cdf,
        ppf=ppf,
    )


def weak_positive_distribution() -> ThrowDistribution:
    """KW kernel (Stage 3 / M2 Gate 1): weak positive mismatch, amplitude ablation.

    Frozen definition (Gate 1; lambda = 1/4 is not a sweep parameter):

        G_W(x) = x + (1/4) x(1-x) = (5/4) x - (1/4) x^2
        g_W(x) = 5/4 - x/2
        Delta_W(x) = G_W(x) - x = (1/4) x(1-x) >= 0

    Same spatial shape and sign as the canonical q_near mismatch
    Delta_C(x) = x(1-x), at 1/4 of its amplitude. The analytic inverse CDF
    follows from solving 5x - x^2 = 4u (small root):
    x = (5 - sqrt(25 - 16u)) / 2.
    """
    def cdf(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        return 1.25 * x - 0.25 * x * x

    def ppf(u: np.ndarray) -> np.ndarray:
        u = np.asarray(u, dtype=float)
        return 0.5 * (5.0 - np.sqrt(25.0 - 16.0 * u))

    return ThrowDistribution(
        name="weak_positive",
        label="G_W(x) = 5/4 x - 1/4 x^2 (weak positive mismatch, lambda = 1/4)",
        pdf=lambda x: 1.25 - 0.5 * np.asarray(x, dtype=float),
        cdf=cdf,
        ppf=ppf,
    )


def flat_top_distribution() -> ThrowDistribution:
    """KS kernel (Stage 3 / M2 Gate 1): positive spatial-shape control.

    Frozen piecewise mismatch (broad flat-top hump instead of the canonical
    parabolic one), matched to canonical q_near in sign, peak, integrated
    mismatch, symmetry and endpoints:

        Delta_S(x) = 3/4 x        , 0 <= x <= 1/3
                     1/4          , 1/3 <= x <= 2/3
                     3/4 (1 - x)  , 2/3 <= x <= 1

        G_S(x) = x + Delta_S(x) = 7/4 x | x + 1/4 | 3/4 + x/4   (same pieces)
        g_S(x) = 7/4 | 1 | 1/4

    Structural properties (verified by tests, exact rational arithmetic):
    Delta_S >= 0, max Delta_S = 1/4 on the flat top, integral_0^1 Delta_S
    = 1/24 + 1/12 + 1/24 = 1/6 — identical to canonical max/integral.
    Continuous at both breakpoints (7/12 and 11/12 respectively).
    """
    def cdf(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        return np.where(
            x <= 1.0 / 3.0,
            1.75 * x,
            np.where(x <= 2.0 / 3.0, x + 0.25, 0.75 + 0.25 * x),
        )

    def ppf(u: np.ndarray) -> np.ndarray:
        u = np.asarray(u, dtype=float)
        return np.where(
            u <= 7.0 / 12.0,
            4.0 * u / 7.0,
            np.where(u <= 11.0 / 12.0, u - 0.25, 4.0 * u - 3.0),
        )

    def pdf(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        return np.where(x <= 1.0 / 3.0, 1.75, np.where(x <= 2.0 / 3.0, 1.0, 0.25))

    return ThrowDistribution(
        name="flat_top",
        label="G_S(x) = x + flat-top mismatch (peak 1/4, integral 1/6)",
        pdf=pdf,
        cdf=cdf,
        ppf=ppf,
    )


def witness_aligned_distribution() -> ThrowDistribution:
    """XA kernel (Stage 3 / M2 Gate 2): Witness-Aligned sign-changing mismatch.

    Frozen piecewise mismatch — positive over most of the canonical
    witness-support region, negative in the Far side (sign boundary 2/3):

        Delta_A(x) = x/2          , 0   <= x <= 1/3
                     1/3 - x/2    , 1/3 <= x <= 2/3
                     2/3 - x      , 2/3 <= x <= 5/6
                     x - 1        , 5/6 <= x <= 1

        G_A(x) = x + Delta_A(x) = 3/2 x | x/2 + 1/3 | 2/3 | 2x - 1
        g_A(x) = 3/2 | 1/2 | 0 | 2   (all >= 0)

    Delta_A > 0 on (0, 2/3), = 0 at x in {0, 2/3, 1}, < 0 on (2/3, 1);
    max Delta_A = +1/6 (x = 1/3), min = -1/6 (x = 5/6). "Witness-Aligned"
    is the experimental ROLE name only, not a claim that rebound must
    occur. XO is the exact pointwise sign reflection: Delta_O = -Delta_A.
    """
    def cdf(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        return np.where(
            x <= 1.0 / 3.0,
            1.5 * x,
            np.where(x <= 2.0 / 3.0, 0.5 * x + 1.0 / 3.0,
                     np.where(x <= 5.0 / 6.0, 2.0 / 3.0 + 0.0 * x,
                              2.0 * x - 1.0)),
        )

    def ppf(u: np.ndarray) -> np.ndarray:
        u = np.asarray(u, dtype=float)
        # left-continuity convention on the flat piece G_A = 2/3 over
        # [2/3, 5/6] (ppf is unused by the deterministic expected-mass path)
        return np.where(
            u <= 0.5,
            2.0 * u / 3.0,
            np.where(u <= 2.0 / 3.0, 2.0 * (u - 1.0 / 3.0),
                     np.where(u < 1.0, (u + 1.0) / 2.0, 1.0 + 0.0 * u)),
        )

    def pdf(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        return np.where(x <= 1.0 / 3.0, 1.5,
                        np.where(x <= 2.0 / 3.0, 0.5,
                                 np.where(x <= 5.0 / 6.0, 0.0, 2.0)))

    return ThrowDistribution(
        name="witness_aligned",
        label="G_A: sign-changing mismatch, Delta_A > 0 on (0, 2/3), "
              "< 0 on (2/3, 1), max |Delta| = 1/6",
        pdf=pdf,
        cdf=cdf,
        ppf=ppf,
    )


def witness_opposed_distribution() -> ThrowDistribution:
    """XO kernel (Stage 3 / M2 Gate 2): Witness-Opposed sign-changing mismatch.

    Exact pointwise sign reflection of XA: Delta_O(x) = -Delta_A(x), so
    abs(Delta_O) == abs(Delta_A) everywhere (the Gate 2 control property);
    G_O = 2x - G_A:

        G_O(x) = x/2            , 0   <= x <= 1/3
                 3/2 x - 1/3    , 1/3 <= x <= 2/3
                 2x - 2/3       , 2/3 <= x <= 5/6
                 1              , 5/6 <= x <= 1

        g_O(x) = 1/2 | 3/2 | 2 | 0   (all >= 0)

    Delta_O < 0 on (0, 2/3) (over the canonical witness-support band),
    > 0 on (2/3, 1) (Far side); max |Delta_O| = 1/6.
    """
    def cdf(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        return np.where(
            x <= 1.0 / 3.0,
            0.5 * x,
            np.where(x <= 2.0 / 3.0, 1.5 * x - 1.0 / 3.0,
                     np.where(x <= 5.0 / 6.0, 2.0 * x - 2.0 / 3.0,
                              1.0 + 0.0 * x)),
        )

    def ppf(u: np.ndarray) -> np.ndarray:
        u = np.asarray(u, dtype=float)
        # left-continuity convention on the flat piece G_O = 1 over [5/6, 1]
        return np.where(
            u <= 1.0 / 6.0,
            2.0 * u,
            np.where(u <= 2.0 / 3.0, 2.0 * (u + 1.0 / 3.0) / 3.0,
                     np.where(u < 1.0, (u + 2.0 / 3.0) / 2.0,
                              5.0 / 6.0 + 0.0 * u)),
        )

    def pdf(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        return np.where(x <= 1.0 / 3.0, 0.5,
                        np.where(x <= 2.0 / 3.0, 1.5,
                                 np.where(x <= 5.0 / 6.0, 2.0, 0.0)))

    return ThrowDistribution(
        name="witness_opposed",
        label="G_O: sign-reflected mismatch, Delta_O = -Delta_A "
              "(negative on the canonical witness band, positive in Far)",
        pdf=pdf,
        cdf=cdf,
        ppf=ppf,
    )


def near_zone_mask(K: int, a: float) -> np.ndarray:
    """Boolean mask of bins fully inside the near zone x < a.

    The deterministic version can only represent the near zone exactly when
    `a` lies on a bin edge; the M0 baseline (a = 0.50, K = 100) satisfies
    this, and anything else is rejected rather than approximated.
    """
    if a <= 0.0 or a >= 1.0:
        raise ValueError(f"near-zone boundary a must lie in (0, 1), got {a}")
    if abs(a * K - round(a * K)) > 1e-9:
        raise ValueError(
            f"a={a} is not on a bin edge for K={K}; the deterministic "
            "version refuses to approximate a straddling bin"
        )
    right_edges = np.arange(1, K + 1) / K
    return right_edges <= a


def uniformity_l2(bin_masses: np.ndarray, K: int) -> float:
    """Primary M0 uniformity error (histogram squared-L2 distance).

        U(t) = sum_{i=1}^{K} ( p_i(t) - 1/K )^2

    where p_i(t) is the fraction of total mass in bin i at round t and the
    ideal fraction of a perfectly uniform state is 1/K. U = 0 exactly for a
    uniform state; smaller is more uniform. The mean squared deviation
    (MSE) reading is U / K.
    """
    if bin_masses.shape != (K,):
        raise ValueError(f"expected {K} bin masses, got shape {bin_masses.shape}")
    dev = bin_masses - 1.0 / K
    return float(np.sum(dev * dev))


def total_variation(bin_masses: np.ndarray, K: int) -> float:
    """Secondary sanity-check metric (total variation distance to uniform).

        TV(t) = 0.5 * sum_{i=1}^{K} | p_i(t) - 1/K |

    Only a sanity check in M0; conclusions must come from the primary U(t).
    """
    return 0.5 * float(np.sum(np.abs(bin_masses - 1.0 / K)))


def noise_floor_l2(N: int, K: int) -> float:
    """Expected U for a perfectly uniform finite sample of N particles.

    For a multinomial histogram of N iid uniform draws,

        E[ sum_i (p_i - 1/K)^2 ] = (K - 1) / (K * N).

    Differences below this level are indistinguishable from finite-sample
    noise in a single Monte Carlo realization.
    """
    return (K - 1) / (K * N)
