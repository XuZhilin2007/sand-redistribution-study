"""M1B.2 micro-scaling diagnostics tests (docs/M1B_MICRO_CORRECTION.md).

Run:  python tests/test_m1b_2_micro_scaling.py

Verifies:
  1. exact action identities from the update rule: moved = alpha*sweep_mass;
     prefix-mass change == -net_export (the net-export definition is exact);
  2. L1_state_step consistency with the state history;
  3. regression anchors: stopping times vs committed M1B results and the
     K=100 trajectory vs committed M1A.1;
  4. the annihilation necessary condition F(a_t) <= 0.5 (q_near, near-end
     excess present) holds at every exact stop;
  5. the headline falsification: late-stage gross action medians are
     K-independent (CV of X far below CV of K*X) in the committed audit CSV;
  6. M1B epsilon stops sit at gross-churn scale (moved*K >> 1), not at
     micro-correction scale.
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
K_LIST = [50, 100, 200, 400]
q = DISTRIBUTIONS["near"]

# --- 1. exact action identities on the canonical dynamics --------------------------
run = run_m1a_deterministic(q, K=K, T_max=500, alpha=ALPHA, stop_tolerance=TOL)
states = np.asarray(run["states"], dtype=float)
ok_moved = ok_export = ok_L1 = True
for t, r in enumerate(run["rows"][:500]):
    if not r["active"]:
        continue
    j_t = int(round(r["a_t"] * K))
    if abs(r["removed_mass"] - ALPHA * r["prefix_mass"]) > 1e-12:
        ok_moved = False
    F_q_a = float(q.cdf(np.array([r["a_t"]]))[0])
    net_export = r["removed_mass"] * (1.0 - F_q_a)
    prefix_before = float(states[t][:j_t].sum())
    prefix_after = float(states[t + 1][:j_t].sum())
    if abs((prefix_after - prefix_before) + net_export) > 1e-10:
        ok_export = False
    L1 = float(np.abs(states[t + 1] - states[t]).sum())
    if L1 < net_export - 1e-12:  # L1 step must dominate the one-way export
        ok_L1 = False
check("moved_mass = alpha * sweep_mass (all sampled rounds)", ok_moved)
check("prefix-mass change == -net_export (exact update identity)", ok_export)
check("L1_state_step >= net_export (export is one component of the step)", ok_L1)

# --- 2. regression anchors ----------------------------------------------------------
m1b_csv = Path(__file__).resolve().parents[1] / "experiments/m1b_tolerance/results/epsilon_k_summary.csv"
m1b_rows = list(csv.DictReader(m1b_csv.open(encoding="utf-8")))
m1b_stop = {(int(r["K"]), float(r["epsilon"])): r["stopping_round"] for r in m1b_rows}
ok = True
for Kk in [50, 100, 200, 400]:
    rk = run_m1a_deterministic(q, K=Kk, T_max=5000, alpha=ALPHA, stop_tolerance=TOL)
    Dm = np.array([r["D_max"] for r in rk["rows"]])
    for eps in [0.0, 0.005, 0.01, 0.02, 0.04]:
        tol_e = TOL if eps == 0.0 else eps
        hit = np.nonzero(Dm <= tol_e)[0]
        t_here = str(int(hit[0])) if len(hit) else "none"
        if t_here != m1b_stop[(Kk, eps)]:
            ok = False
check("stopping times reproduce committed M1B epsilon_k_summary.csv (20 cases)", ok)

m11_csv = Path(__file__).resolve().parents[1] / "experiments/m1a_1_long_horizon/results/long_horizon_trajectory_K100.csv"
m11 = list(csv.DictReader(m11_csv.open(encoding="utf-8")))
run100 = run_m1a_deterministic(q, K=100, T_max=5000, alpha=ALPHA, stop_tolerance=TOL)
same = all(
    abs(run100["rows"][t]["U_L2"] - float(m11[t]["U_L2"]))
    <= 1e-12 * abs(float(m11[t]["U_L2"]))
    and run100["rows"][t]["a_t"] == float(m11[t]["a_t"])
    for t in range(101)
)
check("K=100 trajectory matches committed M1A.1 (t=0..100)", same)

# --- 3. annihilation necessary condition at every exact stop ------------------------
ann_csv = Path(__file__).resolve().parents[1] / "experiments/m1b_2_micro_scaling/results/annihilation_events.csv"
if ann_csv.exists():
    ann_rows = list(csv.DictReader(ann_csv.open(encoding="utf-8")))
    ok_ann = all(r["F_a_le_0.5"] == "True" for r in ann_rows)
    check(f"F(a) <= 0.5 at all {len(ann_rows)} exact stops (necessary condition)",
          ok_ann and len(ann_rows) == 4)
    ok_clear = all(int(r["positive_prefixes_post"]) == 0 for r in ann_rows)
    check("post-stop states have zero positive prefixes", ok_clear)
else:
    check("annihilation CSV present (run the experiment first)", False)

# --- 4. headline falsification from the committed scaling CSV -----------------------
scale_csv = Path(__file__).resolve().parents[1] / "experiments/m1b_2_micro_scaling/results/cross_K_scaling.csv"
if scale_csv.exists():
    srows = list(csv.DictReader(scale_csv.open(encoding="utf-8")))
    for qname in ["moved_mass", "net_export_mass", "L1_state_step"]:
        r = next(x for x in srows if x["quantity"] == qname)
        # honest claim: NO quantity scales as 1/K (|beta| far from 1); X always
        # collapses better than K*X. net_export is noisier (dip-round argmax
        # wanders between small and large prefixes) but beta = -0.216, not -1.
        check(f"{qname}: no 1/K scaling (CV_X < CV_KX, |beta| << 1)",
              float(r["CV_of_X"]) < float(r["CV_of_KX"])
              and abs(float(r["loglog_slope_beta"])) < 0.5,
              f"CV_X={r['CV_of_X']} CV_KX={r['CV_of_KX']} beta={r['loglog_slope_beta']}")
    over_csv = Path(__file__).resolve().parents[1] / "experiments/m1b_2_micro_scaling/results/epsilon_overlay.csv"
    orows = list(csv.DictReader(over_csv.open(encoding="utf-8")))
    ok_regime = all(float(r["moved_mass_times_K"]) > 10 for r in orows)
    check(f"all {len(orows)} M1B stops occur at gross-churn scale (moved*K >> 1)", ok_regime)
else:
    check("scaling CSV present (run the experiment first)", False)

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1B.2 micro-scaling tests passed")
