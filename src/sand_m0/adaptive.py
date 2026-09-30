"""M1A: adaptive sweep boundary — the minimal state-dependent feedback model.

Difference from M0: the sweep region is no longer a fixed x < a. Each round
the model computes the cumulative excess

    D_t(j) = cumulative_mass(bins 1..j) - j/K        (j = 1..K)

which measures how much more sand the prefix [0, j/K] holds than a perfectly
uniform state. The sweep region is the prefix with the largest excess,

    j_t = argmax_j D_t(j),   a_t = j_t / K,   D_max(t) = max_j D_t(j).

If D_max(t) > stop_tolerance, bins 1..j_t are swept (each loses a fraction
alpha of its own mass; the removed mass is resprayed over [0, 1] with the
fixed throw distribution q). If D_max(t) <= stop_tolerance, no prefix holds
a cumulative near-side excess and the model stops sweeping — the model-internal
natural stopping condition (a mathematical rule, NOT a human visual rule).

Because the swept region is chosen from the current state, the update
operator changes from round to round: M0's single-curve description and its
alpha time-rescaling do not carry over (that loss is studied, not assumed).
"""

from __future__ import annotations

import numpy as np

from .model import ThrowDistribution
from .simulate import uniformity_l2, total_variation


def cumulative_excess(mass: np.ndarray, K: int) -> np.ndarray:
    """D(j) for j = 1..K: prefix mass minus the uniform prefix share j/K.

    D(K) = 0 identically (total mass 1), so max_j D(j) >= 0 always; it is
    zero exactly when every prefix holds at most its uniform share.
    """
    if mass.shape != (K,):
        raise ValueError(f"expected {K} bin masses, got shape {mass.shape}")
    return np.cumsum(mass) - np.arange(1, K + 1) / K


def boundary_decision(mass: np.ndarray, K: int) -> tuple[int, float, float]:
    """Adaptive decision: (j_t, a_t, D_max) from the current state.

    Ties in the argmax resolve to the smallest j (deterministic).
    """
    D = cumulative_excess(mass, K)
    j_t = int(np.argmax(D)) + 1
    return j_t, j_t / K, float(D[j_t - 1])


def excess_margin(mass: np.ndarray, K: int) -> float:
    """Gap between the best and second-best prefix excess at the argmax.

    Small margins flag boundary decisions that sampling noise could flip
    (the MC sensitivity motivation); they do not affect the update.
    """
    D = cumulative_excess(mass, K)
    j_t = int(np.argmax(D))
    rest = np.delete(D, j_t)
    return float(D[j_t] - rest.max())


def deficit_fill_redistribution(p_minus: np.ndarray, K: int, removed_mass: float) -> np.ndarray:
    """M1C Variant B kernel: target-aware deficit-filling redistribution.

    Replaces ONLY the respray step of the canonical M1A update. Instead of
    redistributing the removed mass with the fixed throw distribution q, the
    mass is placed back in proportion to each bin's current positive deficit
    against the uniform target u_i = 1/K (the target the cumulative-excess
    selection rule and the uniformity metric already encode):

        deficit_i = max(u_i - p_minus_i, 0)
        add_i     = removed_mass * deficit_i / sum_j deficit_j

    Selection (cumulative-excess argmax), the removal rule (alpha) and the
    stopping rules are untouched.

    Mathematical invariant (exact real arithmetic): sum(p_minus) = 1 -
    removed_mass, so the total deficit equals the post-removal positive
    excess plus removed_mass, hence sum(deficit) >= removed_mass > 0
    whenever a sweep removed mass. Two consequences: the fill can never
    push a deficit bin above u_i (add_i <= deficit_i), and the degenerate
    case sum(deficit) == 0 cannot occur. The explicit errors below turn any
    violation of these identities — i.e. an implementation bug or a
    representation change — into a hard failure instead of a silent
    fallback to the canonical q kernel.
    """
    deficit = np.maximum(1.0 / K - p_minus, 0.0)
    total_deficit = float(deficit.sum())
    if removed_mass <= 0.0:
        return p_minus.copy()
    if total_deficit <= 0.0:
        raise RuntimeError(
            "deficit_fill: total deficit is zero while removed mass > 0 — "
            "violates the mass-conservation identity; refusing to fall back to q"
        )
    if removed_mass > total_deficit * (1.0 + 1e-9):
        raise RuntimeError(
            "deficit_fill: removed mass exceeds total deficit — violates the "
            "mass-conservation identity; refusing to fall back to q"
        )
    return p_minus + removed_mass * (deficit / total_deficit)


def run_m1a_deterministic(
    dist: ThrowDistribution,
    *,
    K: int,
    T_max: int,
    alpha: float,
    stop_tolerance: float,
    redistribution: str = "throw",
    redistribution_dist: ThrowDistribution | None = None,
) -> dict:
    """Deterministic expected-mass M1A trajectory with adaptive boundary.

    The initial state is the throw distribution itself (as in M0). Rows are
    state-indexed (matching run_deterministic): row t carries the metrics of
    state t and the decision for the transition t -> t+1; `active` is True
    iff that sweep is actually executed (never at the horizon t = T_max).
    Once the stopping condition fires the state stays frozen; later rows
    repeat it with active=False. `stopped_at` is None if the excess was
    still positive at every decision up to the horizon.

    `redistribution` selects only how the removed mass re-enters the system
    (M1C ablation); the selection rule, removal strength and stopping rules
    are identical for every value:

      * "throw" (default): canonical M1A/M1B behaviour — removed mass is
        resprayed with the fixed throw distribution q;
      * "deficit_fill": M1C Variant B — removed mass is redistributed by
        `deficit_fill_redistribution` (target-aware deficit filling).

    `redistribution_dist` (M1C.1) optionally splits the two roles the `dist`
    argument plays in the canonical model — initial condition and fixed
    respray law: when given, the initial state still comes from `dist` while
    the "throw" kernel resprays with `redistribution_dist` (Variant U: fixed
    uniform respray, still fixed and deficit-unaware). None (default) keeps
    the canonical single-distribution behaviour bit-identical.
    """
    q = dist.bin_probs(K)
    q_throw = (
        redistribution_dist.bin_probs(K)
        if redistribution_dist is not None
        else q
    )
    mass = q.copy()
    near_zone = np.arange(1, K + 1) / K <= 0.5  # fixed physical zone x < 0.5

    states: list[np.ndarray] = [mass.copy()]
    rows: list[dict] = []
    stopped_at = None
    for t in range(T_max + 1):
        D = cumulative_excess(mass, K)
        j_t = int(np.argmax(D)) + 1  # ties -> smallest index (deterministic)
        a_t = j_t / K
        D_max = float(D[j_t - 1])
        rest = np.delete(D, j_t - 1)
        margin = D_max - float(rest.max())
        would_sweep = stopped_at is None and D_max > stop_tolerance
        active = would_sweep and t < T_max
        prefix_mass = float(mass[:j_t].sum())
        near_mass = float(mass[near_zone].sum())
        removed = alpha * prefix_mass if active else 0.0

        rows.append(
            {
                "t": t,
                "active": active,
                "j_t": j_t,
                "a_t": a_t,
                "D_max": D_max,
                "D_margin": margin,
                "prefix_mass": prefix_mass,
                "near_mass_05": near_mass,
                "removed_mass": removed,
                "U_L2": uniformity_l2(mass, K),
                "TV": total_variation(mass, K),
            }
        )
        if active:
            new_mass = mass.copy()
            new_mass[:j_t] = mass[:j_t] * (1.0 - alpha)
            if redistribution == "throw":
                mass = new_mass + removed * q_throw
            elif redistribution == "deficit_fill":
                mass = deficit_fill_redistribution(new_mass, K, removed)
            else:
                raise ValueError(f"unknown redistribution kernel: {redistribution!r}")
        if not would_sweep and stopped_at is None:
            stopped_at = t
        if t < T_max:
            states.append(mass.copy())

    return {
        "rows": rows,
        "states": states,
        "final_state": states[-1],
        "stopped_at": stopped_at,  # None: still sweeping at T_max
        "redistribution": redistribution,
    }


def run_m1a_mc(
    dist: ThrowDistribution,
    *,
    K: int,
    T_max: int,
    alpha: float,
    stop_tolerance: float,
    N: int,
    seed: int,
) -> dict:
    """Monte Carlo particle M1A: the boundary is chosen from the histogram.

    Spot-check only. The argmax decision is made on a finite-sample histogram,
    so it can flip on sampling noise near ties — recorded, not smoothed.
    """
    rng = np.random.default_rng(seed)
    pos = dist.sample(rng, N)
    states: list[np.ndarray] = []
    rows: list[dict] = []
    stopped_at = None
    for t in range(T_max):
        idx = np.minimum((pos * K).astype(np.int64), K - 1)
        counts = np.bincount(idx, minlength=K)
        p = counts / N
        D = np.cumsum(p) - np.arange(1, K + 1) / K
        j_t = int(np.argmax(D)) + 1
        a_t = j_t / K
        D_max = float(D[j_t - 1])
        active = D_max > stop_tolerance and stopped_at is None
        in_prefix = idx < j_t
        near_mass = float((pos < 0.5).mean())
        prefix_mass = float(in_prefix.mean())
        u = uniformity_l2(p, K)
        if stopped_at is None and D_max <= stop_tolerance:
            stopped_at = t
        if active:
            selected = in_prefix & (rng.random(N) < alpha)
            n_sel = int(selected.sum())
            if n_sel:
                pos[selected] = dist.sample(rng, n_sel)
        rows.append(
            {
                "t": t,
                "active": active,
                "j_t": j_t,
                "a_t": a_t,
                "D_max": D_max,
                "prefix_mass": prefix_mass,
                "near_mass_05": near_mass,
                "removed_frac": (n_sel / N) if active else 0.0,
                "U_L2": u,
                "TV": total_variation(p, K),
            }
        )
    states.append(pos.copy())
    return {"rows": rows, "final_state": states[-1], "stopped_at": stopped_at}
