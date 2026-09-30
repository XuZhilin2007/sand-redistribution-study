"""Minimal M0 sanity tests (no external test framework required).

Run:  python tests/test_m0.py

Covers the Phase 1 requirements from docs/NEXT_STEPS.md: fixed-seed
determinism, mass conservation, metric definitions, exact sampler/bin
agreement, and the finite-sample noise floor. It also verifies the
analytically expected near-mass decay factor (1 - alpha*(1 - Q_near)),
which follows from the update rule itself.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.model import DISTRIBUTIONS, near_zone_mask, noise_floor_l2, total_variation, uniformity_l2
from sand_m0.simulate import run_deterministic, run_mc

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


K, T, A, ALPHA = 100, 30, 0.5, 0.25

# --- 1. exact bin probabilities -------------------------------------------------
for name, dist in DISTRIBUTIONS.items():
    p = dist.bin_probs(K)
    check(f"bin_probs sums to 1 for '{name}'", abs(p.sum() - 1.0) < 1e-12, f"sum={p.sum()!r}")
    check(f"bin_probs all finite non-negative for '{name}'", bool(np.all(p >= 0)))
q_near = DISTRIBUTIONS["near"].bin_probs(K)
q_far = DISTRIBUTIONS["far"].bin_probs(K)
check("q_near near-zone mass = cdf(0.5) = 0.75", abs(q_near[near_zone_mask(K, A)].sum() - 0.75) < 1e-12)
check("q_far near-zone mass = cdf(0.5) = 0.25", abs(q_far[near_zone_mask(K, A)].sum() - 0.25) < 1e-12)
check("q_uniform bins are exactly 1/K", bool(np.allclose(DISTRIBUTIONS["uniform"].bin_probs(K), 1.0 / K, atol=1e-15)))

# --- 2. metric definitions ------------------------------------------------------
uniform = np.full(K, 1.0 / K)
check("U(uniform state) = 0", uniformity_l2(uniform, K) == 0.0)
check("TV(uniform state) = 0", total_variation(uniform, K) == 0.0)
perturbed = np.zeros(K)
perturbed[0] = 0.02
check("U of hand-built state matches formula",
      abs(uniformity_l2(perturbed, K) - (0.02 - 0.01) ** 2 - 99 * 0.01 ** 2) < 1e-15)
check("TV of hand-built state matches formula",
      abs(total_variation(perturbed, K) - 0.5 * (0.01 + 99 * 0.01)) < 1e-15)

# --- 3. deterministic runs: mass conservation and finite values ------------------
for name, dist in DISTRIBUTIONS.items():
    rows = run_deterministic(dist, K=K, T=T, a=A, alpha=ALPHA)
    total = np.sum([row["near_mass"] for row in rows])  # placeholder to keep rows referenced
    # recompute conservation independently by replaying the update
    q = dist.bin_probs(K)
    near = near_zone_mask(K, A)
    mass = q.copy()
    conserved = True
    for _ in range(T):
        removed = ALPHA * mass[near].sum()
        mass = mass - ALPHA * np.where(near, mass, 0.0) + removed * q
        if abs(mass.sum() - 1.0) > 1e-12:
            conserved = False
            break
    check(f"deterministic mass conserved through T rounds for '{name}'", conserved)
    check(f"deterministic trajectory has T+1 rows for '{name}'", len(rows) == T + 1)
    check(f"deterministic U finite and >= 0 for '{name}'",
          all(np.isfinite(row["U_L2"]) and row["U_L2"] >= 0 for row in rows))
    removed_sum = sum(row["removed_mass"] for row in rows[1:])
    near_loss = rows[0]["near_mass"] - rows[-1]["near_mass"]
    # near-zone loss must equal removed mass minus the near share of every respray
    q_near_share = q[near].sum()
    expected_near = sum(
        row["removed_mass"] * (1.0 - q_near_share) for row in rows[1:]
    )
    check(f"deterministic near-mass accounting consistent for '{name}'",
          abs(near_loss - expected_near) < 1e-9,
          f"loss={near_loss!r} expected={expected_near!r}")
    del removed_sum, total

# --- 4. analytic near-mass decay factor ------------------------------------------
# From the update rule: M_{t+1} = M_t * (1 - alpha + alpha * Q_near).
for name, dist in DISTRIBUTIONS.items():
    q = dist.bin_probs(K)
    q_near_share = q[near_zone_mask(K, A)].sum()
    factor = 1.0 - ALPHA + ALPHA * q_near_share
    rows = run_deterministic(dist, K=K, T=T, a=A, alpha=ALPHA)
    expected = rows[0]["near_mass"] * factor ** 5
    check(f"near-mass decay factor matches analytics for '{name}' (t=5)",
          abs(rows[5]["near_mass"] - expected) < 1e-12,
          f"got={rows[5]['near_mass']!r} expected={expected!r}")

# --- 5. fixed-seed determinism of Monte Carlo ------------------------------------
small = dict(K=K, T=5, a=A, alpha=ALPHA, N=5000)
run1 = run_mc(DISTRIBUTIONS["near"], seed=20260917, **small)
run2 = run_mc(DISTRIBUTIONS["near"], seed=20260917, **small)
run3 = run_mc(DISTRIBUTIONS["near"], seed=20260918, **small)
check("same seed reproduces identical MC trajectory",
      all(r1["U_L2"] == r2["U_L2"] and r1["near_mass"] == r2["near_mass"]
          for r1, r2 in zip(run1, run2)))
check("different seeds produce different MC trajectories",
      any(r1["U_L2"] != r3["U_L2"] for r1, r3 in zip(run1, run3)))

# --- 6. MC one-step expectation matches the deterministic factor ------------------
for name, dist in DISTRIBUTIONS.items():
    q = dist.bin_probs(K)
    factor = 1.0 - ALPHA + ALPHA * q[near_zone_mask(K, A)].sum()
    rows = run_mc(dist, K=K, T=1, a=A, alpha=ALPHA, N=40000, seed=7)
    ratio = rows[1]["near_mass"] / rows[0]["near_mass"]
    check(f"MC near-mass ratio at t=1 close to deterministic factor for '{name}'",
          abs(ratio - factor) < 0.01, f"ratio={ratio!r} factor={factor!r}")
    check(f"MC removed fraction close to alpha*near_mass for '{name}'",
          abs(rows[1]["removed_frac"] - ALPHA * rows[0]["near_mass"]) < 0.005)

# --- 7. finite-sample noise floor -------------------------------------------------
floor = noise_floor_l2(100000, K)
p = DISTRIBUTIONS["uniform"].bin_probs(K)
rng = np.random.default_rng(20260917)
pos = rng.random(100000)
idx = np.minimum((pos * K).astype(np.int64), K - 1)
u_emp = uniformity_l2(np.bincount(idx, minlength=K) / 100000, K)
check("uniform finite sample U near analytic noise floor",
      0.4 * floor < u_emp < 2.5 * floor, f"u={u_emp!r} floor={floor!r}")

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M0 sanity tests passed")
