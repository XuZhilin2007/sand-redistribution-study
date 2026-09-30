"""M1C.2 kernel theory: exact one-step discrepancy identities for the
canonical discrete M1 update (docs/M1C2_TARGET_MATCHED_KERNEL_THEORY.md).

Everything here is derived from the canonical update in
`adaptive.run_m1a_deterministic`: bins i = 1..K, mass p >= 0 with sum 1,
uniform target share T(j) = j/K, cumulative excess

    D(j) = F(j) - T(j),   F(j) = sum_{i<=j} p_i,   D(K) = 0,

adaptive boundary j* = argmax_j D(j) (smallest index on ties), and, on an
active round (D_max > stop_tolerance), removal of M = alpha * F(j*) from the
prefix 1..j* (each prefix bin scaled by (1 - alpha)) followed by respraying
the removed mass with a fixed law g (bin probabilities, sum 1).

One active round therefore acts on the cumulative-excess profile as

  prefix part  (j <= j*):  D_new(j) = (1-alpha) D(j)
                                       + alpha [ F(j*) G(j) - T(j) ]
  suffix part  (j >  j*):  D_new(j) = D(j) - M (1 - G(j))

where G(j) = sum_{i<=j} g_i is the respray law's own discrete CDF.
Splitting G against the target isolates the mismatch term

  D_new(j) = D_target_matched_new(j) + M [ G(j) - T(j) ]   (both regions),

i.e. the respray law re-injects, into every prefix discrepancy, exactly
M times the law's OWN cumulative excess D_g(j) = G(j) - T(j) against the
target. For the target-matched law G = T (Variant U) the mismatch vanishes
and the update is strictly D_max-decreasing whenever D_max > 0 (theorem,
proved in the M1C.2 document; verified numerically by
tests/test_m1c2_kernel_theory.py).

Boundary export identity (both regions, any fixed law): the swept prefix's
mass changes by exactly -M (1 - G(j*)).

Continuum notation (D(x) = F(x) - T(x) etc.) is shorthand only; the
canonical model is the discrete one and no continuum limit is established.
"""

from __future__ import annotations

import numpy as np


def prefix_cdf(mass: np.ndarray) -> np.ndarray:
    """F(j) for j = 1..K (cumulative bin masses)."""
    return np.cumsum(mass)


def one_step_excess(
    mass: np.ndarray, j_star: int, alpha: float, respray_bins: np.ndarray
) -> np.ndarray:
    """Exact D-profile of the next state under one active canonical round.

    Applies the region formulas (prefix/suffix) to the CURRENT state; equals
    D(canonically updated mass) up to floating-point evaluation order.
    """
    K = mass.shape[0]
    F = prefix_cdf(mass)
    G = prefix_cdf(respray_bins)
    M = alpha * F[j_star - 1]
    D = F - np.arange(1, K + 1) / K
    T = np.arange(1, K + 1) / K
    D_new = np.empty(K)
    prefix = slice(0, j_star)
    suffix = slice(j_star, K)
    D_new[prefix] = (1.0 - alpha) * D[prefix] + alpha * (F[j_star - 1] * G[prefix] - T[prefix])
    D_new[suffix] = D[suffix] - M * (1.0 - G[suffix])
    return D_new


def target_matched_excess(mass: np.ndarray, j_star: int, alpha: float) -> np.ndarray:
    """D-profile of the next state when the respray law IS the target (G = T).

    For the uniform target this is Variant U's kernel: the mismatch term
    M [G - T] vanishes identically.
    """
    K = mass.shape[0]
    F = prefix_cdf(mass)
    M = alpha * F[j_star - 1]
    D = F - np.arange(1, K + 1) / K
    T = np.arange(1, K + 1) / K
    D_new = np.empty(K)
    prefix = slice(0, j_star)
    suffix = slice(j_star, K)
    D_new[prefix] = (1.0 - alpha) * D[prefix] - alpha * T[prefix] * (1.0 - F[j_star - 1])
    D_new[suffix] = D[suffix] - M * (1.0 - T[suffix])
    return D_new


def mismatch_reinjection(moved_mass: float, respray_bins: np.ndarray) -> np.ndarray:
    """The per-round reinjection profile M [G(j) - T(j)] = M * D_g(j).

    This is the ONLY difference between the actual one-step update and the
    target-matched one, identically in both regions. For q_near against the
    uniform target it is M (j/K)(1 - j/K) >= 0, peaking at M/4 at j = K/2.
    """
    K = respray_bins.shape[0]
    T = np.arange(1, K + 1) / K
    return moved_mass * (prefix_cdf(respray_bins) - T)


def cumulative_mismatch(respray_bins: np.ndarray) -> np.ndarray:
    """Delta(j) = G(j) - T(j) for the uniform target (per-unit-mass profile)."""
    K = respray_bins.shape[0]
    T = np.arange(1, K + 1) / K
    return prefix_cdf(respray_bins) - T


def contraction_margin(mass: np.ndarray, j_star: int, alpha: float) -> np.ndarray:
    """Pointwise target-matched contraction margin gamma(j) = D_max - D_TM'(j).

    gamma(j) is how much positive cumulative excess the target-matched
    one-step map removes at prefix j relative to the CURRENT D_max. By the
    M1C.2 theorem gamma(j) > 0 at every prefix on an active step; on the
    swept prefix (j <= j*) it is bounded below by alpha * D_max.
    """
    D = cumulative_excess_(mass)
    d_max = float(D.max())
    return d_max - target_matched_excess(mass, j_star, alpha)


def rebound_criterion(
    mass: np.ndarray, j_star: int, alpha: float, respray_bins: np.ndarray
) -> dict:
    """Exact one-step D_max increase (rebound) criterion for a fixed law.

    The actual next profile is D'_TM(j) + M Delta(j); the round rebounds iff
    some prefix's reinjection M Delta(j) exceeds its contraction margin
    gamma(j) = D_max - D'_TM(j):

        rebound  <=>  exists j:  M Delta(j) > gamma(j).

    Exact iff in real arithmetic (M1C.3); by the suffix impossibility remark
    any satisfying j must lie in the swept prefix (j <= j*). Returns the
    boolean, the worst reinjection-minus-margin, and the per-prefix vectors.
    """
    gamma = contraction_margin(mass, j_star, alpha)
    M = alpha * float(prefix_cdf(mass)[j_star - 1])
    reinjection = mismatch_reinjection(M, respray_bins)
    excess = reinjection - gamma
    return {
        "rebound": bool((excess > 0.0).any()),
        "worst_excess": float(excess.max()),
        "gamma": gamma,
        "reinjection": reinjection,
        "moved_mass": M,
    }


def cumulative_excess_(mass: np.ndarray) -> np.ndarray:
    """D(j) for the uniform target (local alias to avoid an import cycle)."""
    K = mass.shape[0]
    return prefix_cdf(mass) - np.arange(1, K + 1) / K


def rebound_residual(
    mass: np.ndarray, j_star: int, alpha: float, respray_bins: np.ndarray
) -> np.ndarray:
    """Exact per-prefix rebound residual M*Delta(j) - gamma(j) (M1C.5).

    Direct form (both regions): moved_mass * Delta(j) - gamma(j). On the
    swept prefix (Region 1, j <= j*) this equals the reduced form

        (1-alpha) D(j) - D_max + alpha * h(j/K),

    h(x) = x [F_a (2 - x) - 1],  F_a = F(j*) = a + D_max,

    the identity underlying the M1C.5 boundary-gate theorem (verified
    numerically to <= 1e-15 on the canonical trajectories).
    """
    K = mass.shape[0]
    D = cumulative_excess_(mass)
    d_max = float(D.max())
    F_a = float(mass[:j_star].sum())
    gamma = d_max - target_matched_excess(mass, j_star, alpha)
    moved = alpha * F_a
    return moved * cumulative_mismatch(respray_bins) - gamma


def state_reduced_gate(mass: np.ndarray, j_star: int, alpha: float) -> dict:
    """M1C.5 state-reduced necessary rebound gate (discrete, exact).

    On Region 1 the residual is bounded by alpha (h(j/K) - D_max) with
    h(x) = x [F_a (2 - x) - 1], F_a = F(j*) = a + D_max, because
    (1-alpha) D(j) - D_max <= -alpha D_max. Hence a rebound (residual > 0
    at some prefix) REQUIRES max_j h(j/K) > D_max.

    Closed-form maximization of h on [0, a]:
      * F_a <= 1/2: h <= 0 everywhere -> rebound impossible (exact);
      * F_a > 1/2: interior maximum at x* = 1 - 1/(2 F_a) with value
        (2 F_a - 1)^2 / (4 F_a); x* <= a iff (a + D_max)(1 - a) <= 1/2,
        which holds for every canonical A state (D_max <= 0.25 < sqrt(2)-1).
    The discrete grid maximum max_j h(j/K) differs from the continuum value
    by O(1/K^2) rounding (both returned).
    """
    K = mass.shape[0]
    D = cumulative_excess_(mass)
    d_max = float(D.max())
    F_a = float(mass[:j_star].sum())
    x = np.arange(1, K + 1) / K
    h = x * (F_a * (2.0 - x) - 1.0)
    h_grid_max = float(h[:j_star].max())
    if F_a > 0.5:
        x_star = 1.0 - 1.0 / (2.0 * F_a)
        h_cf = (2.0 * F_a - 1.0) ** 2 / (4.0 * F_a)
    else:
        x_star, h_cf = None, 0.0
    a_crit = (1.0 - d_max + np.sqrt(d_max * (d_max + 2.0))) / 2.0
    return {
        "F_a": F_a, "D_max": d_max, "a": j_star / K,
        "x_star": x_star, "x_star_le_a": (x_star is not None and x_star <= j_star / K),
        "h_grid_max": h_grid_max, "h_closed_form": h_cf,
        "a_crit": float(a_crit),
        "gate_open": j_star / K > float(a_crit),
        "necessary_only": True,
    }


def a_crit_gate(d_max: float) -> float:
    """Analytic boundary gate a_crit(D): rebound impossible for a <= a_crit.

    Closed form (M1C.5): solving (2(a+D)-1)^2 / (4(a+D)) = D for a gives
    a = [(1 - D) +- sqrt(D(D+2))]/2; the large root is the gate (the small
    root lies in the F_a <= 1/2 branch where h <= 0 makes rebound impossible
    anyway). Continuum-style envelope: the discrete grid condition
    max_j h(j/K) > D differs by O(1/K^2) rounding.
    """
    if d_max < 0.0:
        raise ValueError(f"D_max must be >= 0, got {d_max}")
    return float((1.0 - d_max + np.sqrt(d_max * (d_max + 2.0))) / 2.0)
