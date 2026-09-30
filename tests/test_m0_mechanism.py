"""Mechanism tests for M0 (docs/M0_MECHANISM.md).

Run:  python tests/test_m0_mechanism.py

Verifies, against the actual simulator (not just against itself):
  1. exact closed-form state  m(t) = z_t q (near) + lambda_t q (far);
  2. near-mass recurrence M_{t+1} = M_t (1 - alpha(1-Q)) and closed form
     M_t = 0.75 rho^t, including against the committed EXP-M0-A CSV;
  3. exact quadratic U(z) reproduces simulated uniformity error;
  4. the optimum z* (discrete) and the continuous z* = 3/4 candidate,
     including the size of the finite-bin correction;
  5. the time-rescaling property: U(z) does not depend on alpha;
  6. boundary behaviour at alpha = 0 and alpha = 1;
  7. round discretization: which integer round wins and by how much.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.mechanism import (
    analytic_trajectory,
    decay_factor,
    mechanism_constants,
    optimal_round,
    optimal_z,
    quadratic_coefficients,
    u_of_z,
)
from sand_m0.model import DISTRIBUTIONS, uniformity_l2
from sand_m0.simulate import run_deterministic, run_deterministic_states

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


K, T, A_PARAM, BASELINE = 100, 30, 0.5, 0.25
NEAR = DISTRIBUTIONS["near"]
ALPHAS = [0.10, 0.25, 0.50, 0.75, 1.00]
C = mechanism_constants(NEAR, K, A_PARAM)

# --- 1. exact closed-form solution vs the stepped simulator ----------------------
for alpha in [0.0, *ALPHAS]:
    stepped = run_deterministic_states(NEAR, K=K, T=T, a=A_PARAM, alpha=alpha)
    closed = analytic_trajectory(NEAR, K=K, T=T, a=A_PARAM, alpha=alpha)
    err = max(float(np.abs(s - m).max()) for s, m in zip(stepped, closed))
    check(f"closed-form state == stepped state (alpha={alpha})", err < 1e-12,
          f"max err={err:.3e}")

# --- 2. near-mass recurrence ------------------------------------------------------
rows_mem = run_deterministic(NEAR, K=K, T=T, a=A_PARAM, alpha=BASELINE)
rho = decay_factor(BASELINE, C.Q_near)
check("rho = 15/16 exactly", rho == 15 / 16, f"rho={rho!r}")
err_mem = max(abs(r["near_mass"] - 0.75 * rho**r["t"]) for r in rows_mem)
check("M_t = 0.75*(15/16)^t vs in-memory simulator", err_mem < 1e-12,
      f"max err={err_mem:.3e}")

csv_path = Path(__file__).resolve().parents[1] / "experiments/m0_baseline/results/EXP-M0-A_deterministic.csv"
rows_csv = list(csv.DictReader(csv_path.open(encoding="utf-8")))
err_csv = max(abs(float(r["near_mass"]) - 0.75 * rho**int(r["t"])) for r in rows_csv)
check("M_t closed form vs committed EXP-M0-A CSV", err_csv <= 5e-11,
      f"max err={err_csv:.3e} (CSV prints 10 decimals, so ~5e-11 is formatting)")

# per-round recurrence identity on in-memory rows
err_rec = max(
    abs(rows_mem[t + 1]["near_mass"] - rows_mem[t]["near_mass"] * (1 - BASELINE * (1 - C.Q_near)))
    for t in range(T)
)
check("recurrence M_{t+1} = M_t*(1 - alpha*(1-Q)) round by round", err_rec < 1e-12,
      f"max err={err_rec:.3e}")

# --- 3. exact quadratic U(z) ------------------------------------------------------
A_q, B_q, C_q = quadratic_coefficients(C)
err_u = max(abs(u_of_z(rho**r["t"], C) - r["U_L2"]) for r in rows_mem)
check("U(z) quadratic reproduces simulated U(t) (baseline)", err_u < 1e-12,
      f"max err={err_u:.3e}")
err_u_csv = max(
    abs(u_of_z(rho**int(r["t"]), C) - float(r["U_L2"])) for r in rows_csv
)
check("U(z) quadratic reproduces committed CSV U column", err_u_csv < 1e-12,
      f"max err={err_u_csv:.3e}")

# identity U(z) = S_n z^2 + S_f ((1-Qz)/Q_f)^2 - 1/K evaluated at arbitrary z
for z in [0.9, 0.772476196, 0.6, 0.3]:
    m = np.where(C.near_mask, z * C.q, ((1 - C.Q_near * z) / C.Q_far) * C.q)
    direct = uniformity_l2(m, K)
    check(f"quadratic U({z:.4f}) matches direct histogram evaluation",
          abs(u_of_z(z, C) - direct) < 1e-14, f"diff={abs(u_of_z(z, C) - direct):.3e}")

# --- 4. optimum -------------------------------------------------------------------
z_star = optimal_z(C)
check("discrete z* within 1e-4 of the continuous 3/4",
      abs(z_star - 0.75) < 1e-4, f"z*={z_star!r} (correction {z_star - 0.75:+.3e})")
u_star = u_of_z(z_star, C)
u_traj = [u_of_z(rho**t, C) for t in range(T + 1)]
check("no trajectory round beats the continuous optimum U(z*)",
      all(u >= u_star - 1e-15 for u in u_traj),
      f"min over rounds={min(u_traj):.6e} vs U(z*)={u_star:.6e}")
t_star = optimal_round(rho, z_star)
check("continuous t* near 4.46 (baseline)", abs(t_star - 4.46) < 0.05, f"t*={t_star:.4f}")

# --- 5. time-rescaling: U(z) independent of alpha ---------------------------------
max_coeff_diff = 0.0
for alpha in ALPHAS:
    c_a = quadratic_coefficients(mechanism_constants(NEAR, K, A_PARAM))
    max_coeff_diff = max(max_coeff_diff, abs(c_a[0] - A_q), abs(c_a[1] - B_q), abs(c_a[2] - C_q))
check("quadratic coefficients identical across all alpha", max_coeff_diff == 0.0,
      f"max coeff diff={max_coeff_diff:.3e}")

# pooled (z, U) points from all simulated alphas lie on one curve
pool_err = 0.0
for alpha in ALPHAS:
    r_a = run_deterministic(NEAR, K=K, T=T, a=A_PARAM, alpha=alpha)
    rho_a = decay_factor(alpha, C.Q_near)
    pool_err = max(pool_err, max(abs(u_of_z(rho_a**r["t"], C) - r["U_L2"]) for r in r_a))
check("pooled (z, U) from all alphas on a single exact curve", pool_err < 1e-12,
      f"max residual={pool_err:.3e}")

# --- 6. boundaries ----------------------------------------------------------------
rows0 = run_deterministic(NEAR, K=K, T=T, a=A_PARAM, alpha=0.0)
check("alpha=0: state never moves (U(t)=U(0) for all t)",
      all(r["U_L2"] == rows0[0]["U_L2"] for r in rows0)
      and all(r["near_mass"] == rows0[0]["near_mass"] for r in rows0))
rows1 = run_deterministic(NEAR, K=K, T=T, a=A_PARAM, alpha=1.0)
check("alpha=1: z_1 = 0.75 exactly (near mass 0.5625 at t=1)",
      abs(rows1[1]["near_mass"] - 0.75 * 0.75) < 1e-12,
      f"near_mass(1)={rows1[1]['near_mass']!r}")

# --- 7. round discretization ------------------------------------------------------
z4, z5 = rho**4, rho**5
check("t=4 is closer to z* than t=5 (hence U(4) < U(5))",
      abs(z4 - z_star) < abs(z5 - z_star) and u_of_z(z4, C) < u_of_z(z5, C),
      f"|z4-z*|={abs(z4 - z_star):.6f} |z5-z*|={abs(z5 - z_star):.6f}")
gap = u_of_z(z4, C) - u_star
check("discretization gap U(4) - U(z*) is positive and small (< 1e-4)",
      0 < gap < 1e-4, f"gap={gap:.3e}")

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M0 mechanism tests passed")
