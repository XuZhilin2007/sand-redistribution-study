"""M1B.1 metric-scale audit tests (docs/M1B_METRIC_SCALE_AUDIT.md).

Run:  python tests/test_m1b_1_metric_scale.py

Verifies:
  1. coarse aggregation identities (mass preserved; B=K reduces exactly to
     the fine-grid U_density);
  2. stopping times reproduce the committed M1B epsilon_k_summary.csv;
  3. K=100 reference trajectory matches the committed M1A.1 trajectory;
  4. retention bounds and argmin correctness on real audit data;
  5. retention math on a hand-built example.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.adaptive import run_m1a_deterministic
from sand_m0.model import DISTRIBUTIONS

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


K, T_MAX, ALPHA, TOL = 100, 5000, 0.25, 1e-12
B_LIST = [5, 10, 25, 50]
K_LIST = [50, 100, 200, 400]
EPSILONS = [0.0, 0.005, 0.01, 0.02, 0.04]

# --- 1. aggregation identity on a synthetic state --------------------------------
m = np.random.default_rng(7).dirichlet(np.ones(100))
m = m / m.sum()
for B in [5, 10, 25, 50, 100]:
    P = m.reshape(B, K // B).sum(axis=1)
    check(f"B={B}: coarse masses sum to 1", abs(P.sum() - 1.0) < 1e-12)
    u_coarse = B * ((P - 1.0 / B) ** 2).sum()
    if B == K:
        u_fine = K * ((m - 1.0 / K) ** 2).sum()
        check("B=K reduces exactly to the fine-grid U_density", u_coarse == u_fine)

# --- 2. committed M1B stopping times reproduced ------------------------------------
m1b_csv = Path(__file__).resolve().parents[1] / "experiments/m1b_tolerance/results/epsilon_k_summary.csv"
m1b_rows = list(csv.DictReader(m1b_csv.open(encoding="utf-8")))
m1b_stop = {(int(r["K"]), float(r["epsilon"])): r["stopping_round"] for r in m1b_rows}
ok = True
for K in K_LIST:
    run = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=T_MAX, alpha=ALPHA,
                                stop_tolerance=TOL)
    Dm = np.array([r["D_max"] for r in run["rows"]])
    for eps in EPSILONS:
        tol_e = TOL if eps == 0.0 else eps
        hit = np.nonzero(Dm <= tol_e)[0]
        t_here = str(int(hit[0])) if len(hit) else "none"
        if t_here != m1b_stop[(K, eps)]:
            ok = False
check("stopping times reproduce committed M1B epsilon_k_summary.csv (20 cases)", ok)

# --- 3. K=100 trajectory matches committed M1A.1 -----------------------------------
m11_csv = Path(__file__).resolve().parents[1] / "experiments/m1a_1_long_horizon/results/long_horizon_trajectory_K100.csv"
m11 = list(csv.DictReader(m11_csv.open(encoding="utf-8")))
run100 = run_m1a_deterministic(DISTRIBUTIONS["near"], K=100, T_max=T_MAX, alpha=ALPHA,
                               stop_tolerance=TOL)  # K explicit: loop above shadows K
ok = all(
    abs(run100["rows"][t]["U_L2"] - float(m11[t]["U_L2"])) <= 1e-12 * abs(float(m11[t]["U_L2"]))
    and run100["rows"][t]["a_t"] == float(m11[t]["a_t"])
    for t in range(101)
)
check("K=100 reference trajectory matches committed M1A.1 (t=0..100)", ok)

# --- 4. audit-data invariants (from the committed audit CSV) ------------------------
audit_csv = Path(__file__).resolve().parents[1] / "experiments/m1b_1_metric_scale/results/stop_quality_by_scale.csv"
rows = list(csv.DictReader(audit_csv.open(encoding="utf-8")))
ok_bounds = all(
    -1e-9 <= float(r["improvement_retention"]) <= 1.0 + 1e-9
    and float(r["improvement_loss"]) >= -1e-9
    and float(r["absolute_gap"]) >= -1e-15
    for r in rows
)
check(f"all {len(rows)} (K,eps,B) stop-quality rows: retention in [0,1], loss >= 0",
      ok_bounds and len(rows) == 95)  # 20 rows at K=50 (B=50==K deduped) + 25 x 3
ok_argmin = all(
    float(r["U_stop"]) >= float(r["U_best"]) - 1e-15 for r in rows
)
check("U_stop >= U_best for every (K,eps,B) (argmin correctness)", ok_argmin)

# --- 5. retention math on a hand-built example --------------------------------------
u_init, u_best, u_stop = 0.3333, 0.003821, 0.008540
loss = (u_stop - u_best) / (u_init - u_best)
retention = 1.0 - loss
check("hand-built example: retention ~0.9857 (matches the small-denominator story)",
      abs(retention - 0.9857) < 1e-3, f"retention={retention:.6f}")

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1B.1 metric-scale tests passed")
