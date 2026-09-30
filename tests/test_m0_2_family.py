"""Phase 2A tests: general-q mechanism and the frozen power family.

Run:  python tests/test_m0_2_family.py

Verifies (all against the canonical simulator, not just against itself):
  1. general-q shape preservation / one-variable trajectory / exact
     quadratic U(z) — including a q OUTSIDE the frozen family;
  2. power family sampling and bin probabilities;
  3. the interior-optimum condition  Q_f*S_n > Q*S_f  (c_n > c_f) and the
     predicted behaviour in all three regimes;
  4. family closed forms  c_n/c_f = (2^(2p+1)-1)/(2^(p+1)-1)  and
     z*_cont = 2(1-u)/(3-4u),  u = 2^-(p+1);
  5. predicted best rounds for the six frozen cases;
  6. time rescaling for general q (alpha only reparameterizes time);
  7. regression anchor: p = 1 reproduces the M0/M0.1 baseline numbers.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.mechanism import (
    decay_factor,
    mechanism_constants,
    optimal_round,
    optimal_z,
    quadratic_coefficients,
    u_of_z,
)
from sand_m0.model import DISTRIBUTIONS, power_distribution, uniformity_l2
from sand_m0.simulate import run_deterministic, run_deterministic_states

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


K, T, A_PARAM, ALPHA0 = 100, 30, 0.5, 0.25
FROZEN_P = [-0.25, 0.0, 0.25, 0.5, 1.0, 2.0]

# --- 1. general q outside the family: shape lock + exact quadratic ---------------
def piecewise_cdf(x):
    x = np.asarray(x, dtype=float)
    tot = 1.15
    out = np.where(
        x < 0.5,
        (1.6 * x - 0.6 * x * x) / tot,
        (0.65 + (0.4 * x + 0.4 * x * x - 0.3)) / tot,
    )
    return np.clip(out, 0.0, 1.0)


class PiecewiseDist:
    """Minimal stand-in distribution outside the power family."""

    name = "piecewise"
    label = "piecewise test q"

    @staticmethod
    def bin_probs(k):
        edges = np.linspace(0.0, 1.0, k + 1)
        p = np.diff(piecewise_cdf(edges))
        return p / p.sum()

    @staticmethod
    def sample(rng, n):
        u = rng.random(n)
        lo = np.zeros(n)
        hi = np.ones(n)
        for _ in range(80):
            mid = (lo + hi) / 2
            lo = np.where(piecewise_cdf(mid) < u, mid, lo)
            hi = np.where(piecewise_cdf(mid) < u, hi, mid)
        return (lo + hi) / 2


C_pw = mechanism_constants(PiecewiseDist(), K, A_PARAM)  # type: ignore[arg-type]
rho_pw = decay_factor(ALPHA0, C_pw.Q_near)
sts = run_deterministic_states(PiecewiseDist(), K=K, T=T, a=A_PARAM, alpha=ALPHA0)  # type: ignore[arg-type]
err_state, err_u = 0.0, 0.0
for t, m_sim in enumerate(sts):
    z = rho_pw**t
    lam = (1.0 - C_pw.Q_near * z) / C_pw.Q_far
    m_th = np.where(C_pw.near_mask, z * C_pw.q, lam * C_pw.q)
    err_state = max(err_state, float(np.abs(m_sim - m_th).max()))
    err_u = max(err_u, abs(uniformity_l2(m_sim, K) - u_of_z(z, C_pw)))
check("general q (outside family): closed-form state matches simulator",
      err_state < 1e-12, f"max err={err_state:.3e}")
check("general q (outside family): U(z) quadratic matches simulator",
      err_u < 1e-12, f"max err={err_u:.3e}")

# --- 2. power family sampling and bin probabilities -------------------------------
for p in FROZEN_P:
    dist = power_distribution(p)
    probs = dist.bin_probs(K)
    check(f"p={p}: bin probs sum to 1 and are non-negative",
          abs(probs.sum() - 1.0) < 1e-12 and bool(np.all(probs >= 0)))
    u = np.array([0.1, 0.37, 0.5, 0.83, 0.999])
    x = dist.ppf(u)
    check(f"p={p}: ppf inverts cdf",
          bool(np.all(np.abs(dist.cdf(x) - u) < 1e-10)))
C = {p: mechanism_constants(power_distribution(p), K, A_PARAM) for p in FROZEN_P}
check("p=0 is exactly the uniform distribution",
      bool(np.allclose(C[0.0].q, 1.0 / K, atol=1e-15)))
check("p=1 is exactly the M0 baseline q_near",
      bool(np.allclose(C[1.0].q, DISTRIBUTIONS["near"].bin_probs(K), atol=1e-15)))

# --- 3. interior condition and the three regimes ----------------------------------
for p in FROZEN_P:
    Cc = C[p]
    Q, Qf, Sn, Sf = Cc.Q_near, Cc.Q_far, Cc.S_near, Cc.S_far
    z_star = optimal_z(Cc)
    # tolerance handles the exact boundary (uniform): fp residues ~1e-19 there
    cond_holds = (Qf * Sn - Q * Sf) > 1e-12
    is_interior = z_star < 1.0 - 1e-9
    check(f"p={p}: interior condition Q_f*S_n > Q*S_f <=> z* < 1",
          cond_holds == is_interior,
          f"cond={cond_holds}, z*={z_star:.5f}")
    rows = run_deterministic(power_distribution(p), K=K, T=T, a=A_PARAM, alpha=ALPHA0)
    best_t = min(range(T + 1), key=lambda t: rows[t]["U_L2"])
    if is_interior:
        improves = best_t > 0
        worsens = rows[-1]["U_L2"] > rows[best_t]["U_L2"]
        check(f"p={p}: interior regime observed (improves then worsens)",
              improves and worsens, f"best_t={best_t}")
    else:
        check(f"p={p}: boundary/counterexample regime observed (t=0 optimal, monotone)",
              best_t == 0 and all(
                  rows[t + 1]["U_L2"] >= rows[t]["U_L2"] - 1e-15 for t in range(T)
              ), f"best_t={best_t}")
check("boundary case p=0: z* = 1 and c_n/c_f = 1 to fp precision",
      abs(optimal_z(C[0.0]) - 1.0) < 1e-9
      and abs((C[0.0].S_near / C[0.0].Q_near) / (C[0.0].S_far / C[0.0].Q_far) - 1.0) < 1e-9)

# --- 4. family closed forms --------------------------------------------------------
for p in FROZEN_P:
    Cc = C[p]
    ratio_disc = (Cc.S_near / Cc.Q_near) / (Cc.S_far / Cc.Q_far)
    ratio_cf = (2.0 ** (2 * p + 1) - 1.0) / (2.0 ** (p + 1) - 1.0)
    # The closed form is the continuous (a=0.5) prediction; the discrete 100-bin
    # implementation carries finite-bin corrections that grow with within-bin
    # concentration (largest on the far-concentrated side, e.g. ~1.6% at p=-0.25).
    tol = 0.03 if p < 0 else 5e-3
    check(f"p={p}: discrete c_n/c_f matches continuous closed form within finite-bin correction",
          abs(ratio_disc - ratio_cf) / ratio_cf < tol,
          f"disc={ratio_disc:.5f} cf={ratio_cf:.5f}")
    u = 2.0 ** (-(p + 1.0))
    z_cont = 2.0 * (1.0 - u) / (3.0 - 4.0 * u)
    z_star = optimal_z(Cc)
    check(f"p={p}: discrete z* close to closed form 2(1-u)/(3-4u)",
          abs(z_star - z_cont) < 0.02, f"disc={z_star:.5f} cont={z_cont:.5f}")

# --- 5. predicted best rounds for the frozen cases ---------------------------------
expected_best = {-0.25: 0, 0.0: 0, 0.25: 1, 0.5: 2, 1.0: 4, 2.0: 11}
for p in FROZEN_P:
    Cc = C[p]
    rho = decay_factor(ALPHA0, Cc.Q_near)
    A, B, Cq = quadratic_coefficients(Cc)
    us = [A * (rho**t) ** 2 + B * (rho**t) + Cq for t in range(T + 1)]
    pred_t = min(range(T + 1), key=lambda t: us[t])
    rows = run_deterministic(power_distribution(p), K=K, T=T, a=A_PARAM, alpha=ALPHA0)
    obs_t = min(range(T + 1), key=lambda t: rows[t]["U_L2"])
    check(f"p={p}: predicted best round {pred_t} == observed {obs_t} == frozen {expected_best[p]}",
          pred_t == obs_t == expected_best[p], f"pred={pred_t} obs={obs_t}")

# --- 6. time rescaling for general q -----------------------------------------------
for p in [0.5, 2.0, -0.25]:
    Cc = C[p]
    A, B, Cq = quadratic_coefficients(Cc)
    worst = 0.0
    for alpha in [0.1, 0.5, 1.0]:
        rows = run_deterministic(power_distribution(p), K=K, T=T, a=A_PARAM, alpha=alpha)
        rho = decay_factor(alpha, Cc.Q_near)
        worst = max(worst, max(abs(r["U_L2"] - (A * (rho**r["t"]) ** 2 + B * (rho**r["t"]) + Cq))
                               for r in rows))
    check(f"p={p}: all spot-check alphas on one exact quadratic U(z)", worst < 1e-12,
          f"max residual={worst:.3e}")

# alpha=1 for the strong-bias case: z_1 = Q = 0.875
rows1 = run_deterministic(power_distribution(2.0), K=K, T=2, a=A_PARAM, alpha=1.0)
check("p=2, alpha=1: near mass after one round = Q*z = 0.875*0.875",
      abs(rows1[1]["near_mass"] - 0.875**2) < 1e-12, f"{rows1[1]['near_mass']!r}")

# --- 7. regression anchor p=1 -------------------------------------------------------
check("p=1 reproduces M0.1 discrete z* = 0.74997187",
      abs(optimal_z(C[1.0]) - 0.7499718732420776) < 1e-12,
      f"z*={optimal_z(C[1.0])!r}")
rho1 = decay_factor(ALPHA0, C[1.0].Q_near)
check("p=1 reproduces M0.1 continuous t* = 4.4581",
      abs(optimal_round(rho1, optimal_z(C[1.0])) - 4.4581) < 1e-3)
check("p=1 near-mass factor is exactly 15/16", rho1 == 15 / 16)

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M0.2 family tests passed")
