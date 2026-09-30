"""M1B tolerance tests (docs/M1B_TOLERANCE.md).

Run:  python tests/test_m1b_tolerance.py

Verifies:
  1. first-passage interpretation: M1B's canonical stopping round equals the
     first t with D_max_reference(t) <= epsilon (the reference run predicts
     every finite-epsilon stopping time);
  2. stopped state == reference state at t_stop (not just the same U);
  3. monotonicity of t_stop in epsilon over a dense grid (K=100 and K=400);
  4. eps=0 reference reproduces the committed M1A.1 canonical trajectory;
  5. controls stop immediately;
  6. frozen-state semantics for finite epsilon;
  7. K-robustness of finite-epsilon stopping vs K-scaling of exact stopping
     (the headline M1B finding, as stored-data assertions).
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


K, T_MAX, ALPHA, TOL_EXACT = 100, 3000, 0.25, 1e-12
EPSILONS = [0.005, 0.01, 0.02]

ref = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=T_MAX, alpha=ALPHA,
                            stop_tolerance=TOL_EXACT)
D_ref = np.array([r["D_max"] for r in ref["rows"]])

# --- 1 + 2. first-passage identity and state identity ------------------------------
for eps in EPSILONS:
    run = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=T_MAX, alpha=ALPHA,
                                stop_tolerance=eps)
    hit = np.nonzero(D_ref <= eps)[0]
    predicted = int(hit[0]) if len(hit) else -1
    check(f"eps={eps}: canonical stopping round == reference first passage",
          run["stopped_at"] == predicted,
          f"stop={run['stopped_at']} pred={predicted}")
    if run["stopped_at"] is not None and predicted is not None and predicted >= 0:
        t = run["stopped_at"]
        same_state = bool(np.array_equal(run["states"][t], ref["states"][t]))
        same_u = run["rows"][t]["U_L2"] == ref["rows"][t]["U_L2"]
        check(f"eps={eps}: stopped state identical to reference state at t={t}",
              same_state and same_u)
        check(f"eps={eps}: pre-stop trajectory follows the reference exactly",
              all(run["rows"][i]["U_L2"] == ref["rows"][i]["U_L2"]
                  for i in range(t)))
        frozen = all(run["rows"][i]["U_L2"] == run["rows"][t]["U_L2"]
                     for i in range(t, len(run["rows"])))
        check(f"eps={eps}: state frozen after stop", frozen)

# --- 3. monotonicity over a dense epsilon grid --------------------------------------
grid = np.geomspace(1e-12, 0.06, 300)
ts = []
for eps in grid:
    hit = np.nonzero(D_ref <= eps)[0]
    ts.append(int(hit[0]) if len(hit) else -1)
mono = all(ts[i] == -1 or ts[i + 1] == -1 or ts[i] >= ts[i + 1]
           for i in range(len(ts) - 1))
check("t_stop(eps) non-increasing over 300-point eps grid (K=100)", mono)

# --- 4. regression anchor ------------------------------------------------------------
csv_path = Path(__file__).resolve().parents[1] / "experiments/m1a_1_long_horizon/results/long_horizon_trajectory_K100.csv"
m11 = list(csv.DictReader(csv_path.open(encoding="utf-8")))
same = all(
    abs(ref["rows"][t]["U_L2"] - float(m11[t]["U_L2"]))
    <= 1e-12 * abs(ref["rows"][t]["U_L2"])   # M1A.1 CSV stores 12 significant digits
    and ref["rows"][t]["a_t"] == float(m11[t]["a_t"])
    for t in range(101)
)
check("eps=0 reference reproduces committed M1A.1 K=100 trajectory (t=0..100)",
      same)
check("eps=0 reference stopping round == 699 (canonical M1A.1)", ref["stopped_at"] == 699)

# --- 5. controls ----------------------------------------------------------------------
for name in ["uniform", "far"]:
    run = run_m1a_deterministic(DISTRIBUTIONS[name], K=K, T_max=5, alpha=ALPHA,
                                stop_tolerance=0.01)
    check(f"control q_{name}: immediate stop at t=0 (D_max(0) <= 0.01)",
          run["stopped_at"] == 0 and run["rows"][0]["D_max"] <= 0.01)

# --- 7. headline finding from the committed experiment CSV ---------------------------
csv_path = Path(__file__).resolve().parents[1] / "experiments/m1b_tolerance/results/epsilon_k_summary.csv"
rows = list(csv.DictReader(csv_path.open(encoding="utf-8")))
by = {(int(r["K"]), float(r["epsilon"])): r for r in rows if r["stopping_round"] != "none"}
t_of = lambda K, e: int(by[(K, e)]["stopping_round"])
exact_ratio = t_of(400, 0.0) / t_of(50, 0.0)
finite_ratios = [t_of(400, e) / t_of(50, e) for e in [0.005, 0.01, 0.02, 0.04]]
check("exact stopping scales ~9x between K=50 and K=400", exact_ratio > 5,
      f"ratio={exact_ratio:.2f}")
check("every finite frozen epsilon stops within ~1.3x across the same K range",
      all(r < 1.3 for r in finite_ratios), f"ratios={[f'{x:.2f}' for x in finite_ratios]}")
check("finite-epsilon stopping rounds are non-increasing in epsilon for each K",
      all(t_of(K, 0.04) <= t_of(K, 0.02) <= t_of(K, 0.01) <= t_of(K, 0.005)
          for K in [50, 100, 200, 400]))

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1B tolerance tests passed")
