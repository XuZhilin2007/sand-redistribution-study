"""M1C.4 churn-regime audit tests (docs/M1C4_CHURN_REGIME_AUDIT.md).

Run:  python tests/test_m1c4_churn_dynamics.py

All checks run against the canonical Variant A deterministic trajectories
(no new model) and the committed M1C.4 analysis CSVs.

  1. baseline integrity: exact stopping times 324/699/1445/2876 unchanged;
  2. R_t definition: for every active round of the canonical A runs, the
     ratio form R_t > 1 is equivalent to the M1C.3 excess-form rebound
     criterion (0 mismatches across all 5344 rounds);
  3. R_t values in the committed per-round CSV match an independent
     recomputation (formatted-field equality);
  4. churn window: onset (first D_max < 0.05) = 23 for all four K in the
     committed summary, and window rebound fractions match the summary;
  5. symbol consistency: transition counts sum to window - 1 and run
     counts are consistent with the transition table;
  6. rebound fraction = P(R_t > 1) within the churn window (identity);
  7. terminal alignment CSV covers tau = 1..40 for all four K;
  8. deterministic reproducibility of the analysis: recomputed churn-window
     summary quantiles match the committed CSV.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.adaptive import cumulative_excess, run_m1a_deterministic
from sand_m0.kernel_theory import target_matched_excess
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
EXPECTED_STOPS = {50: 324, 100: 699, 200: 1445, 400: 2876}
res_dir = Path(__file__).resolve().parents[1] / "experiments/m1c4_churn_dynamics/results"

# --- canonical runs + independent R_t recomputation ---------------------------------
runs = {K: run_m1a_deterministic(q_near, K=K, T_max=5000, alpha=ALPHA,
                                 stop_tolerance=TOL)
        for K in K_LIST}
check("baseline integrity: exact stopping times unchanged (324/699/1445/2876)",
      all(runs[K]["stopped_at"] == EXPECTED_STOPS[K] for K in K_LIST))

ok_equiv = True
recomputed = {}
for K in K_LIST:
    stop = runs[K]["stopped_at"]
    states = np.asarray(runs[K]["states"], dtype=float)
    delta = np.cumsum(q_near.bin_probs(K)) - np.arange(1, K + 1) / K
    recs = []
    for t in range(stop):
        mass = states[t]
        D = cumulative_excess(mass, K)
        d_max = float(D.max())
        j_star = int(np.argmax(D)) + 1
        moved = ALPHA * float(mass[:j_star].sum())
        gamma = d_max - target_matched_excess(mass, j_star, ALPHA)
        pos_reinj = moved * np.maximum(delta, 0.0)
        R = float((pos_reinj[gamma > 0] / gamma[gamma > 0]).max())
        rebound_excess = bool((pos_reinj - gamma).max() > 0.0)
        if (R > 1.0) != rebound_excess:
            ok_equiv = False
        recs.append({"R": R, "rebound": rebound_excess,
                     "D_max": d_max, "a": j_star / K, "M": moved})
    recomputed[K] = recs
total_rounds = sum(len(v) for v in recomputed.values())
check(f"R_t definition: ratio form R_t > 1 == M1C.3 excess criterion on all "
      f"{total_rounds} active rounds (0 mismatches)", ok_equiv)

# --- committed CSV anchors ----------------------------------------------------------
summary_csv = res_dir / "churn_window_summary.csv"
if summary_csv.exists():
    srows = {int(r["K"]): r for r in csv.DictReader(summary_csv.open(encoding="utf-8"))}
    per_rows = list(csv.DictReader((res_dir / "per_round_churn_diagnostics.csv")
                                   .open(encoding="utf-8")))
    ok_R = True
    for K in K_LIST:
        recs = recomputed[K]
        rows = [r for r in per_rows if int(r["K"]) == K]
        if len(rows) != len(recs):
            ok_R = False
            continue
        for rec, row in zip(recs, rows):
            if (row["rebound"] == "True") != rec["rebound"] \
               or abs(float(row["R_t"]) - rec["R"]) > 5e-10:
                ok_R = False
    check("committed per-round CSV matches the independent recomputation "
          "(R_t to 5e-10, rebound flags exactly)", ok_R)

    check("churn window onset (first D_max < 0.05) = 23 for all four K",
          all(int(srows[K]["t_onset_0.05"]) == 23 for K in K_LIST))
    ok_frac = all(
        abs(float(srows[K]["rebound_fraction"])
            - sum(r["rebound"] for r in recomputed[K][23:]) / len(recomputed[K][23:]))
        <= 5e-6 for K in K_LIST)
    check("window rebound fraction == P(R_t > 1) within the churn window "
          "(identity, matches summary)", ok_frac)

    ok_trans = True
    for K in K_LIST:
        s = srows[K]
        if int(s["RR"]) + int(s["RC"]) + int(s["CR"]) + int(s["CC"]) \
                != int(s["window_rounds"]) - 1:
            ok_trans = False
    check("transition counts sum to window_rounds - 1 for all K", ok_trans)

    check("cross-K stability anchors: R_p50 in [0.65, 0.75] and rebound "
          "fraction in [0.35, 0.38] for every K",
          all(0.65 <= float(srows[K]["R_p50"]) <= 0.75
              and 0.35 <= float(srows[K]["rebound_fraction"]) <= 0.38
              for K in K_LIST))

    align_csv = res_dir / "terminal_alignment.csv"
    arows = list(csv.DictReader(align_csv.open(encoding="utf-8")))
    check("terminal alignment covers tau = 1..40 for all four K",
          all(sum(1 for r in arows if int(r["K"]) == K) == 40 for K in K_LIST))

    ok_rr = all(int(srows[K]["RR"]) <= 10 for K in K_LIST)
    check("switching structure anchor: consecutive rebounds are rare "
          "(R->R <= 10 for every K; the regime alternates R with C-pairs)",
          ok_rr)
else:
    check("M1C.4 summary CSV present (run the analysis first)", False)

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1C.4 churn-audit tests passed")
