"""M0 simulation drivers: deterministic expected-mass and Monte Carlo runs."""

from __future__ import annotations

import numpy as np

from .model import (
    ThrowDistribution,
    near_zone_mask,
    total_variation,
    uniformity_l2,
)


def run_deterministic_states(
    dist: ThrowDistribution, *, K: int, T: int, a: float, alpha: float
) -> list[np.ndarray]:
    """Bin-mass state vector at every round t = 0..T.

    Same dynamics as run_deterministic, but returns the full K-bin mass
    vector of each round for mechanism analyses (exact-solution checks,
    rescaling studies). The state sequence is bit-identical to the one
    run_deterministic steps through internally.
    """
    q = dist.bin_probs(K)
    near = near_zone_mask(K, a)
    mass = q.copy()
    states = [mass.copy()]
    for _ in range(T):
        removed = alpha * float(mass[near].sum())
        near_part = np.where(near, mass, 0.0)
        mass = mass - alpha * near_part + removed * q
        states.append(mass.copy())
    return states


def run_deterministic(
    dist: ThrowDistribution, *, K: int, T: int, a: float, alpha: float
) -> list[dict]:
    """Expected-mass (infinite-particle) M0 trajectory.

    Tracks the exact expected mass of each of K bins. The initial state is
    the throw distribution itself. Each round the near-zone mass shrinks by
    a factor alpha and the removed mass is redistributed over the bin
    probabilities of the same throw distribution. Total mass stays 1.

    Returns one row per round t = 0..T with the state at t and the mass
    removed during the transition (t-1) -> t (0.0 at t = 0).
    """
    near = near_zone_mask(K, a)
    states = run_deterministic_states(dist, K=K, T=T, a=a, alpha=alpha)

    rows: list[dict] = []
    for t, mass in enumerate(states):
        near_mass = float(mass[near].sum())
        removed_mass = 0.0 if t == 0 else alpha * float(states[t - 1][near].sum())
        rows.append(
            {
                "t": t,
                "near_mass": near_mass,
                "removed_mass": removed_mass,
                "U_L2": uniformity_l2(mass, K),
                "TV": total_variation(mass, K),
            }
        )
    return rows


def run_mc(
    dist: ThrowDistribution,
    *,
    K: int,
    T: int,
    a: float,
    alpha: float,
    N: int,
    seed: int,
) -> list[dict]:
    """Monte Carlo particle M0 trajectory for one seed.

    N particles carry positions in [0, 1]. Each round, particles with
    x < a are selected independently with probability alpha and each
    selected particle is resampled independently from the throw
    distribution. Metrics are computed from the K-bin histogram of
    positions; the near-zone mass uses exact positions (x < a).
    """
    rng = np.random.default_rng(seed)
    pos = dist.sample(rng, N)

    rows: list[dict] = []
    removed_frac = 0.0
    for t in range(T + 1):
        idx = np.minimum((pos * K).astype(np.int64), K - 1)
        p = np.bincount(idx, minlength=K) / N
        rows.append(
            {
                "seed": seed,
                "t": t,
                "near_mass": float((pos < a).mean()),
                "removed_frac": removed_frac,
                "U_L2": uniformity_l2(p, K),
                "TV": total_variation(p, K),
            }
        )
        if t < T:
            selected = (pos < a) & (rng.random(N) < alpha)
            n_selected = int(selected.sum())
            removed_frac = n_selected / N
            if n_selected:
                pos[selected] = dist.sample(rng, n_selected)
    return rows


def best_round(trajectories: list[dict], key: str = "U_L2") -> tuple[int, float]:
    """Round with the minimal metric value, given full trajectories.

    Ties resolve to the earliest round. Works for a single trajectory
    (list of per-round rows) and ignores the extra columns.
    """
    ts = [row["t"] for row in trajectories]
    values = [row[key] for row in trajectories]
    best_index = min(range(len(values)), key=lambda i: values[i])
    return ts[best_index], values[best_index]
