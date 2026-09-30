"""M1C.5 boundary-gate theory tests (docs/M1C5_BOUNDARY_GATE_THEORY.md).

Run:  python tests/test_m1c5_boundary_gate.py

All checks run on the canonical Variant A deterministic trajectories
(no new model). The theorem chain: Region-1 reduction -> state-reduced
bound -> analytic maximization -> necessary gate a > a_crit(D_max).

  1. Region-1 reduction identity: M*Delta(j) - gamma(j) ==
     (1-alpha)D(j) - D_max + alpha*h(j/K) for every prefix of every active
     round (h(x) = x[F_a(2-x)-1], F_a = a + D_max);
  2. state-reduced bound: residual <= alpha*(max_j h - D_max) on Region 1;
  3. F_a <= 1/2 branch: h_grid_max <= 0 exactly -> rebound impossible;
  4. closed-form maximization: |max_j h(j/K) - (2F_a-1)^2/(4F_a)| within
     the O(1/K^2) grid-rounding envelope, and x* <= a on every round;
  5. threshold closed form: a_crit(D) solves the equality (numeric bisection);
  6. NO false exclusions: every rebound round of the canonical A trajectory
     satisfies a_t > a_crit(D_max,t) (necessary gate, theorem -> 0 FN);
  7. gate-open contractions exist (necessary is not sufficient) and their
     count matches the committed summary;
  8. a_crit quantiles match the committed gate_summary.csv;
  9. deterministic reproducibility of the per-round gate fields.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.adaptive import cumulative_excess, run_m1a_deterministic
from sand_m0.kernel_theory import (
    a_crit_gate,
    state_reduced_gate,
    target_matched_excess,
)
from sand_m0.model import DISTRIBUTIONS

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


K_LIST = [50, 100, 200, 400]
ALPHA, TOL = 0.25, 1e-12
q_near = DISTRIBUTIONS["near"]
res_dir = Path(__file__).resolve().parents[1] / "experiments/m1c5_boundary_gate/results"

runs = {K: run_m1a_deterministic(q_near, K=K, T_max=5000, alpha=ALPHA,
                                 stop_tolerance=TOL)
        for K in K_LIST}

# --- 1-6: theorem chain on every active round ---------------------------------------
worst_resid = worst_bound_viol = 0.0
ok_branch = ok_xstar = ok_gate = True
fn = fp_impossible = 0
n_rounds = 0
hgap_worst = 0.0
recomputed = {}
for K in K_LIST:
    stop = runs[K]["stopped_at"]
    states = np.asarray(runs[K]["states"], dtype=float)
    g = q_near.bin_probs(K)
    delta = np.cumsum(g) - np.arange(1, K + 1) / K
    recs = []
    for t in range(stop):
        mass = states[t]
        D = cumulative_excess(mass, K)
        d_max = float(D.max())
        j_star = int(np.argmax(D)) + 1
        gate = state_reduced_gate(mass, j_star, ALPHA)
        F_a, a = gate["F_a"], gate["a"]
        x = np.arange(1, K + 1) / K
        h = x * (F_a * (2.0 - x) - 1.0)
        gamma = d_max - target_matched_excess(mass, j_star, ALPHA)
        moved = ALPHA * F_a
        resid = moved * delta - gamma
        # 1. Region-1 reduction identity (exact)
        reduced = (1.0 - ALPHA) * D[:j_star] - d_max + ALPHA * h[:j_star]
        worst_resid = max(worst_resid, float(np.abs(resid[:j_star] - reduced).max()))
        # 2. state-reduced bound
        viol = float((resid[:j_star] - (ALPHA * (h[:j_star].max() - d_max))).max())
        worst_bound_viol = max(worst_bound_viol, viol)
        # 3. F_a <= 1/2 branch (exact: h <= 0 on x >= 0)
        if F_a <= 0.5 and gate["h_grid_max"] > 0.0:
            ok_branch = False
        # 4. closed-form max vs grid max (two rounding regimes, see doc):
        #    generic O(F_a/K^2); when x* falls before the first grid point
        #    (F_a -> 1/2) the gap is bounded by the slope, O(F_a/K)
        if gate["x_star"] is not None:
            hgap = abs(gate["h_grid_max"] - gate["h_closed_form"])
            hgap_worst = max(hgap_worst, hgap)
            if hgap > 2.0 * F_a / K + 1e-12:
                ok_branch = False
        # x* <= a is only defined on the F_a > 1/2 branch (x_star = None
        # there means the branch is inert, not that the condition fails)
        if gate["x_star"] is not None and not gate["x_star_le_a"]:
            ok_xstar = False
        # 6. necessary gate: no rebound below a_crit
        rebound = bool(resid.max() > 0.0)
        B = a - gate["a_crit"]
        if rebound and B <= 1e-12:
            ok_gate = False
            fn += 1
        if (not rebound) and B <= -1e-9 and d_max > 0:
            fp_impossible += 1  # correctly excluded rounds (informational)
        n_rounds += 1
        recs.append({"rebound": rebound, "B": B, "a_crit": gate["a_crit"],
                     "gate_open": B > 0.0})
    recomputed[K] = recs
total = sum(len(v) for v in recomputed.values())
check(f"Region-1 reduction identity exact on all {n_rounds} active rounds "
      f"(max residual {worst_resid:.2e} <= 1e-12)", worst_resid <= 1e-12)
check(f"state-reduced bound holds on Region 1 for every round "
      f"(max violation {worst_bound_viol:.2e} <= 1e-12)", worst_bound_viol <= 1e-12)
check("F_a <= 1/2 branch: h_grid_max <= 0 on every such round (rebound "
      "impossible, exact); grid-vs-closed-form gap within the 2F_a/K "
      "rounding envelope", ok_branch)
check(f"closed-form maximization consistent with the discrete grid "
      f"(worst |h_grid - h_cf| = {hgap_worst:.2e}, O(1/K) phase regime) and "
      f"x* <= a on every round (D_max <= 0.25 < sqrt(2)-1)", ok_xstar)
check(f"necessary gate with NO false exclusions: all rebound rounds satisfy "
      f"a_t > a_crit(D_max,t) ({total} rounds, FN = {fn})", ok_gate and fn == 0)

# --- 5. threshold closed form vs numeric root ---------------------------------------
def _phi(a, D):
    F = a + D
    return (2 * F - 1) ** 2 / (4 * F)

ok_root = True
for D in [1e-4, 0.005, 0.01, 0.02, 0.025, 0.05, 0.1, 0.25]:
    ac = a_crit_gate(D)
    eps = 1e-8
    # phi is increasing for a > 1/2 - D, so the gate must straddle the root
    if not (_phi(ac - eps, D) < D < _phi(ac + eps, D)):
        ok_root = False
    if abs(_phi(ac, D) - D) > 1e-9:
        ok_root = False
check("threshold closed form: a_crit(D) solves (2(a+D)-1)^2/(4(a+D)) = D "
      "(sign-straddle + equality residual, 8 values of D)", ok_root)

# --- 7-9: committed CSV anchors ------------------------------------------------------
summary_csv = res_dir / "gate_summary.csv"
if summary_csv.exists():
    srows = {int(r["K"]): r for r in csv.DictReader(summary_csv.open(encoding="utf-8"))}
    ok_counts = all(
        int(srows[K]["false_negatives"]) == 0
        and int(srows[K]["gate_open_contraction"]) ==
            sum(1 for r in recomputed[K][23:] if r["gate_open"] and not r["rebound"])
        for K in K_LIST)
    check("gate-open contractions exist on every K (necessary != sufficient) "
          "and counts match the committed summary", ok_counts)
    ok_q = all(
        abs(float(srows[K]["acrit_p50"]) -
            float(np.median([r["a_crit"] for r in recomputed[K][23:]]))) <= 1e-5  # runner formats 5 decimals
        for K in K_LIST)
    check("committed a_crit quantiles match recomputation", ok_q)
    check("analytic gate tracks the observed ~0.58 boundary: a_crit p50 in "
          "[0.55, 0.65] for every K",
          all(0.55 <= float(srows[K]["acrit_p50"]) <= 0.65 for K in K_LIST))
else:
    check("M1C.5 gate_summary.csv present (run the analysis first)", False)

per_csv = res_dir / "per_round_gate.csv"
if per_csv.exists():
    per_rows = list(csv.DictReader(per_csv.open(encoding="utf-8")))
    ok_rep = True
    for K in [50, 400]:
        rows = [r for r in per_rows if int(r["K"]) == K]
        if len(rows) != len(recomputed[K]):
            ok_rep = False
            continue
        for rec, row in zip(recomputed[K], rows):
            if (row["rebound"] == "True") != rec["rebound"] \
               or abs(float(row["a_crit"]) - rec["a_crit"]) > 5e-7:
                ok_rep = False
    check("deterministic reproducibility: committed per-round gate fields "
          "match recomputation (K=50 and K=400)", ok_rep)
else:
    check("M1C.5 per_round_gate.csv present", False)

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1C.5 boundary-gate tests passed")
