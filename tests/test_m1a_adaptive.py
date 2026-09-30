"""M1A adaptive-boundary tests (docs/M1A_MODEL.md).

Run:  python tests/test_m1a_adaptive.py

Verifies:
  1. cumulative excess definition on hand-built states (incl. D(K)=0);
  2. first M1A decision for the three control distributions (uniform and
     far stop immediately; q_near starts at a=0.5 exactly);
  3. the first M1A update coincides with the M0 update (a_0 = 0.5);
  4. closed-form second boundary  x* = (2L-1)/(2L),  L = (1-Qz)/Q_f  for
     every spot-check alpha, and its monotonicity in alpha;
  5. mass conservation and non-negativity along adaptive trajectories;
  6. stopping semantics: frozen state, no sweep after stop, stopped_at
     recorded once;
  7. alpha time-rescaling breaks: a_1 (second decision) differs across
     alphas -> genuinely different trajectories;
  8. fixed-seed determinism of the MC variant.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.adaptive import cumulative_excess, run_m1a_deterministic, run_m1a_mc
from sand_m0.model import DISTRIBUTIONS
from sand_m0.simulate import run_deterministic

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


K, T_MAX, ALPHA0, TOL = 100, 100, 0.25, 1e-12
ALPHAS = [0.10, 0.25, 0.50, 1.00]

# --- 1. cumulative excess definition ----------------------------------------------
m = np.full(K, 1.0 / K)
check("uniform state: all D(j)=0", bool(np.allclose(cumulative_excess(m, K), 0.0, atol=1e-15)))
m2 = np.zeros(K)
m2[0] = 0.02
D2 = cumulative_excess(m2, K)
check("hand-built state: D(1)=0.02-0.01, then D decreases by 1/K per empty bin",
      abs(D2[0] - 0.01) < 1e-15 and abs(D2[1] - D2[0] + 1.0 / K) < 1e-15)
m3 = np.zeros(K)
m3[-1] = 1.0
check("all mass at far end: every prefix under-filled (D<=0)",
      bool(np.all(cumulative_excess(m3, K) <= 1e-15)) and
      abs(cumulative_excess(m3, K)[-1]) < 1e-15)

# --- 2. first decision per distribution --------------------------------------------
for name, expect_stop in [("uniform", True), ("far", True), ("near", False)]:
    res = run_m1a_deterministic(DISTRIBUTIONS[name], K=K, T_max=3, alpha=ALPHA0,
                                stop_tolerance=TOL)
    first = res["rows"][0]
    if expect_stop:
        check(f"q_{name}: stops immediately at t=0 (D_max<=tol, no sweep)",
              not first["active"] and res["stopped_at"] == 0 and first["D_max"] <= TOL,
              f"D_max={first['D_max']:.2e}")
    else:
        check("q_near: first decision a=0.5 exactly (D=x(1-x) peak)",
              first["active"] and first["a_t"] == 0.5 and
              abs(first["D_max"] - 0.25) < 1e-12,
              f"a={first['a_t']} D_max={first['D_max']}")

# --- 3. first M1A update == first M0 update -----------------------------------------
m1a = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=2, alpha=ALPHA0,
                            stop_tolerance=TOL)
m0_rows = run_deterministic(DISTRIBUTIONS["near"], K=K, T=1, a=0.5, alpha=ALPHA0)
check("M1A state at t=1 identical to M0 state at t=1 (same first update)",
      abs(m1a["rows"][1]["U_L2"] - m0_rows[1]["U_L2"]) < 1e-15
      and abs(m1a["rows"][1]["near_mass_05"] - m0_rows[1]["near_mass"]) < 1e-15,
      f"U: {m1a['rows'][1]['U_L2']!r} vs {m0_rows[1]['U_L2']!r}")

# --- 4. closed-form second boundary and alpha monotonicity --------------------------
expected_a1 = {0.10: 0.53, 0.25: 0.58, 0.50: 0.64, 1.00: 0.71}
prev = 0.0
mono = True
for alpha in ALPHAS:
    res = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=3, alpha=alpha,
                                stop_tolerance=TOL)
    a1 = res["rows"][1]["a_t"]
    z = 1.0 - alpha * 0.25
    L = (1.0 - 0.75 * z) / 0.25
    x_star = (2.0 * L - 1.0) / (2.0 * L)
    check(f"alpha={alpha}: second decision a_1={a1} matches closed form "
          f"x*=({2*L:.4f}-1)/(2*{L:.4f})={x_star:.5f}",
          abs(a1 - expected_a1[alpha]) < 1e-12 and abs(a1 - x_star) < 0.01,
          f"a_1={a1}")
    mono = mono and a1 > prev
    prev = a1
check("second decision boundary monotone increasing in alpha", mono)

# --- 5. conservation and non-negativity ---------------------------------------------
res = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=T_MAX, alpha=ALPHA0,
                            stop_tolerance=TOL)
ok_mass = all(abs(s.sum() - 1.0) < 1e-12 for s in res["states"])
ok_nonneg = all(bool(np.all(s >= -1e-15)) for s in res["states"])
check(f"M1A q_near: mass conserved through {T_MAX} rounds", ok_mass)
check("M1A q_near: bin masses non-negative throughout", ok_nonneg)

# --- 6. stopping semantics -----------------------------------------------------------
res_u = run_m1a_deterministic(DISTRIBUTIONS["uniform"], K=K, T_max=5, alpha=ALPHA0,
                              stop_tolerance=TOL)
frozen = all(
    res_u["rows"][t]["U_L2"] == res_u["rows"][4]["U_L2"] for t in range(5)
) and all(not r["active"] for r in res_u["rows"])
check("uniform control: state frozen after t=0 stop, no sweep rows",
      frozen and res_u["stopped_at"] == 0)

# --- 7. alpha time-rescaling breaks ---------------------------------------------------
a_seqs = {}
u_at_2 = {}
for alpha in ALPHAS:
    r = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=10, alpha=alpha,
                              stop_tolerance=TOL)
    a_seqs[alpha] = [row["a_t"] for row in r["rows"]]
    u_at_2[alpha] = r["rows"][2]["U_L2"]
first_div = {}
for alpha in ALPHAS:
    for t in range(T_MAX // 10):
        if abs(a_seqs[alpha][t] - a_seqs[0.25][t]) > 1e-12:
            first_div[alpha] = t
            break
check("boundary sequences diverge at the SECOND decision (t=1) for every other alpha",
      all(v == 1 for v in first_div.values()) and len(first_div) == len(ALPHAS) - 1,
      f"first divergence: {first_div}")
check("U at t=2 genuinely differs across alphas (not a time-rescaled copy)",
      len({round(v, 15) for v in u_at_2.values()}) == len(ALPHAS),
      f"U(t=2): {u_at_2}")

# --- 8. MC determinism ----------------------------------------------------------------
r1 = run_m1a_mc(DISTRIBUTIONS["near"], K=K, T_max=5, alpha=ALPHA0,
                stop_tolerance=TOL, N=20000, seed=20260917)
r2 = run_m1a_mc(DISTRIBUTIONS["near"], K=K, T_max=5, alpha=ALPHA0,
                stop_tolerance=TOL, N=20000, seed=20260917)
check("MC M1A: same seed reproduces identical trajectory",
      all(a["U_L2"] == b["U_L2"] and a["a_t"] == b["a_t"]
          for a, b in zip(r1["rows"], r2["rows"])))

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1A adaptive tests passed")
