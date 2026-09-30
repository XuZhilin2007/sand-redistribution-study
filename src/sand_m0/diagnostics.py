"""Diagnostics for M1A.1: long-horizon behaviour, boundary ambiguity, recurrence.

These tools are analysis-only: they never change the M1A rule, the sweep
policy or any state. Everything here reads states produced by
`adaptive.run_m1a_deterministic` and computes:

  * the exact identity between the discrete cumulative excess maximum and
    the continuous one-sided CDF discrepancy (one-sided Kolmogorov–Smirnov
    functional against the uniform target);
  * boundary-decision ambiguity diagnostics (value gap, spatial separation
    of the top-two prefixes, near-optimal boundary sets);
  * L1 recurrence statistics over lags (limit-cycle / near-repetition probe);
  * rolling-window summary statistics.
"""

from __future__ import annotations

import numpy as np


def sup_cdf_minus_uniform(mass: np.ndarray, K: int, fine: int = 200001) -> float:
    """sup_{x in [0,1]} (F(x) - x) for the step-function CDF of `mass`.

    For a right-continuous step CDF with jumps at the bin edges j/K, the
    supremum of F(x) - x is attained at the left endpoints of the bins,
    i.e. exactly at the grid points j/K. Hence this continuous supremum is
    EQUAL to max_j D(j) (the discrete cumulative-excess maximum), which the
    function verifies numerically on a fine grid. F(x) - x against the
    uniform target CDF (x) is the one-sided Kolmogorov–Smirnov functional.
    """
    edges = np.arange(K + 1) / K
    F = np.concatenate([[0.0], np.cumsum(mass)])
    grid = np.linspace(0.0, 1.0, fine)
    F_grid = F[np.searchsorted(edges, grid, side="right") - 1]
    return float(np.max(F_grid - grid))


def boundary_diagnostics(mass: np.ndarray, K: int) -> dict:
    """Boundary-decision ambiguity diagnostics for one state.

    Returns the top-two prefixes by cumulative excess D(j) (the winner
    excluded when searching the runner-up), their value gap and spatial
    separation, and the near-optimal boundary sets
    A_p = {j : D(j) >= p * D_max} (p = 0.99, 0.95) with their spatial span,
    number of connected components and the width of the component holding
    the argmax. Diagnostics only; the sweep policy always uses the argmax.
    """
    D = np.cumsum(mass) - np.arange(1, K + 1) / K
    j1 = int(np.argmax(D))
    v1 = float(D[j1])
    D_wo = D.copy()
    D_wo[j1] = -np.inf
    j2 = int(np.argmax(D_wo))
    v2 = float(D_wo[j2])

    out: dict = {
        "j1": j1 + 1,
        "a1": (j1 + 1) / K,
        "D_max": v1,
        "j2": j2 + 1,
        "a2": (j2 + 1) / K,
        "D_second": v2,
        "value_gap": v1 - v2,
        "spatial_separation": abs(j1 - j2) / K,
    }
    for label, p in [("1pct", 0.99), ("5pct", 0.95)]:
        # A_p = {j : D(j) >= p * D_max}. Meaningful when D_max > 0 (the only
        # regime where M1A sweeps); for D_max <= 0 the set degenerates.
        if v1 <= 0:
            out[f"span_{label}"] = 0.0
            out[f"components_{label}"] = 1
            out[f"main_width_{label}"] = 1.0 / K
            out[f"farthest_point_{label}"] = 0.0
            continue
        mask = D >= p * v1
        idx = np.nonzero(mask)[0]
        # connected components = maximal runs of consecutive indices
        breaks = np.nonzero(np.diff(idx) > 1)[0]
        n_comp = int(len(breaks) + 1)
        comp_id = np.zeros(len(idx), dtype=int)
        comp_id[breaks + 1] = np.arange(1, n_comp)
        pos = int(np.searchsorted(idx, j1))
        comp = comp_id[pos]
        members = idx[comp_id == comp]
        out[f"span_{label}"] = float((idx[-1] - idx[0] + 1) / K)
        out[f"components_{label}"] = n_comp
        out[f"main_width_{label}"] = float((members[-1] - members[0] + 1) / K)
        # farthest near-optimal boundary from the argmax: near 0 = the
        # near-optimal set hugs the winner; large = genuinely far options
        out[f"farthest_point_{label}"] = float(np.abs(idx - j1).max()) / K
    return out


def window_stats(values: np.ndarray) -> dict:
    """min / max / mean / std (population) of a 1-D window of values."""
    v = np.asarray(values, dtype=float)
    return {
        "min": float(v.min()),
        "max": float(v.max()),
        "mean": float(v.mean()),
        "std": float(v.std()),
    }


def recurrence_l1(states: list[np.ndarray], max_lag: int, second_half: bool = True) -> dict:
    """L1 recurrence probe: for each lag k, min and mean of ||m_t - m_{t-k}||_1.

    Only rounds t in the second half of the trajectory are used (a settled
    regime, if any, must live there). Returns {k: {"min", "mean"}}. Small
    minima indicate near-repeated states (limit-cycle-like behaviour);
    uniformly large minima indicate no short-period recurrence. This probe
    does NOT establish or refute chaos.
    """
    n = len(states)
    start = n // 2 if second_half else 0
    arr = np.asarray(states[start:], dtype=float)  # (m, K)
    out: dict = {}
    for k in range(1, max_lag + 1):
        if arr.shape[0] <= k:
            break
        d = np.abs(arr[k:] - arr[:-k]).sum(axis=1)
        out[k] = {"min": float(d.min()), "mean": float(d.mean())}
    return out


def largest_jumps(a_seq: list[float], n: int = 10) -> list[dict]:
    """The n largest round-to-round boundary jumps |a_t - a_{t-1}|."""
    jumps = [
        {"t": t, "a_before": a_seq[t - 1], "a_after": a_seq[t],
         "jump": abs(a_seq[t] - a_seq[t - 1])}
        for t in range(1, len(a_seq))
    ]
    jumps.sort(key=lambda r: r["jump"], reverse=True)
    return jumps[:n]
