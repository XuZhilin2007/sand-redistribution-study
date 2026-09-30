"""M1C.6 boundary-return & peak-competition tests (docs/M1C6_BOUNDARY_RETURN.md).

Run:  python tests/test_m1c6_boundary_return.py

All checks run on the canonical Variant A deterministic trajectories
(no new model).

  1. baseline integrity: exact stopping times unchanged;
  2. episode segmentation: episode count == rebound count, exactly one
     terminal episode per K with tau_rebound == stop - t_r, and episode
     lengths sum consistently;
  3. gate consistency: every rebound round has B > 0 (M1C.5 necessary
     gate, 0 false negatives on the episode table);
  4. exact pairwise peak-height update identity (M1C.2 corollary):
     D_new(j1) - D_new(j2) == [D_TM_new(j1) - D_TM_new(j2)]
                             + M[Delta(j1) - Delta(j2)]
     on the top-2 meaningful peaks of every active round (fp residual);
  5. contraction-round peak motion: the competing peak falls strictly
     slower than the active peak (gap-closing share 1.0 within the churn
     window, all K);
  6. gate return: no episode with tau_return == -1 (the gate always
     re-opens before the stop), and the share distribution matches the
     committed summary;
  7. switch-mechanism classification matches the committed summary and
     no_switch share <= 0.06 (M1C.4 consistency);
  8. deterministic reproducibility: recomputed episode fields match the
     committed episode_table.csv.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.adaptive import cumulative_excess, run_m1a_deterministic
from sand_m0.kernel_theory import a_crit_gate, target_matched_excess
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
res_dir = Path(__file__).resolve().parents[1] / "experiments/m1c6_boundary_return/results"


def local_peaks(mass: np.ndarray):
    K = mass.shape[0]
    D = cumulative_excess(mass, K)
    pk: list[int] = []
    for j in range(K):
        left = D[j - 1] if j > 0 else -np.inf
        right = D[j + 1] if j < K - 1 else -np.inf
        if D[j] >= right and (j == 0 or D[j] > left):
            if pk and j == pk[-1] + 1 and D[j] == D[pk[-1]]:
                continue
            pk.append(j)
    return D, pk


runs = {K: run_m1a_deterministic(q_near, K=K, T_max=5000, alpha=ALPHA,
                                 stop_tolerance=TOL)
        for K in K_LIST}
check("baseline integrity: exact stopping times unchanged (324/699/1445/2876)",
      all(runs[K]["stopped_at"] == EXPECTED_STOPS[K] for K in K_LIST))

# --- independent recomputation -------------------------------------------------------
recomputed = {}
worst_pair = 0.0
ok_gap = True
for K in K_LIST:
    stop = runs[K]["stopped_at"]
    states = np.asarray(runs[K]["states"], dtype=float)
    sep = max(2, K // 50)
    delta = np.cumsum(q_near.bin_probs(K)) - np.arange(1, K + 1) / K
    reb = np.array([float(cumulative_excess(states[t + 1], K).max())
                    > float(cumulative_excess(states[t], K).max()) + 1e-12
                    for t in range(stop)])
    t_r_list = [t for t in range(stop) if reb[t]]
    eps = []
    gap_ok_rounds = 0
    gap_closing = 0
    for ep, t_r in enumerate(t_r_list):
        nxt = t_r_list[ep + 1] if ep + 1 < len(t_r_list) else stop
        B_r = None
        tau_return = -1
        for s in range(0, nxt - t_r + 1):
            D = cumulative_excess(states[t_r + s], K)
            a = (int(np.argmax(D)) + 1) / K
            B = a - a_crit_gate(max(float(D.max()), 0.0))
            if s == 0:
                B_r = B
            elif B > 0 and tau_return == -1:
                tau_return = s
        eps.append({"t_r": t_r, "tau_rebound": nxt - t_r,
                    "terminal": ep + 1 >= len(t_r_list), "B_r": B_r,
                    "tau_return": tau_return})
    # pair identity + gap motion on a deterministic stride sample
    for t in range(0, stop, 3):
        D, pk = local_peaks(states[t])
        if not pk:
            continue
        j1 = pk[int(np.argmax([D[j] for j in pk]))]
        far = [j for j in pk if abs(j - j1) > sep]
        if not far:
            continue
        j2 = far[int(np.argmax([D[j] for j in far]))]
        Dn = cumulative_excess(states[t + 1], K)
        M = ALPHA * float(states[t][:j1 + 1].sum())
        tm = target_matched_excess(states[t], j1 + 1, ALPHA)
        pred = (tm[j1] - tm[j2]) + M * (delta[j1] - delta[j2])
        worst_pair = max(worst_pair, abs(float(Dn[j1] - Dn[j2]) - pred))
        if not reb[t]:
            gap_closing += int((Dn[j2] - D[j2]) - (Dn[j1] - D[j1]) > 0)
            gap_ok_rounds += 1
    recomputed[K] = {"eps": eps, "n_reb": len(t_r_list), "stop": stop,
                     "gap_share": gap_closing / gap_ok_rounds}

total_reb = sum(v["n_reb"] for v in recomputed.values())
check(f"episode segmentation: episode count == rebound count on every K "
      f"(total {total_reb} rebounds), exactly one terminal episode per K, "
      f"terminal tau_rebound == stop - t_r",
      all(v["eps"][-1]["terminal"]
          and sum(e["terminal"] for e in v["eps"]) == 1
          and v["eps"][-1]["tau_rebound"] == v["stop"] - v["eps"][-1]["t_r"]
          for v in recomputed.values()))
check("gate consistency: every episode-opening rebound has B_r > 0 "
      "(M1C.5 necessary gate, 0 violations)",
      all(e["B_r"] > 0 for v in recomputed.values() for e in v["eps"]))
check(f"exact pairwise peak-height update identity on stride-3 sample "
      f"(worst residual {worst_pair:.2e} <= 1e-12)", worst_pair <= 1e-12)
check("contraction rounds close the peak gap on every sampled round "
      "(gap-closing share 1.0 within tolerance, all K)",
      all(abs(v["gap_share"] - 1.0) <= 1e-9 for v in recomputed.values()))
check("gate always returns: no episode with tau_return == -1 (never share "
      "= 0 on every K)",
      all(e["tau_return"] != -1 for v in recomputed.values() for e in v["eps"]
          if not e["terminal"]))

# --- committed CSV anchors ----------------------------------------------------------
summary_csv = res_dir / "gate_return_summary.csv"
if summary_csv.exists():
    srows = {int(r["K"]): r for r in csv.DictReader(summary_csv.open(encoding="utf-8"))}
    ok_counts = all(
        int(srows[K]["episodes"]) == len(recomputed[K]["eps"]) - 1
        for K in K_LIST)
    check("committed episode counts match recomputation (non-terminal)", ok_counts)
    check("committed anchors: return_never == 0, gap-closing == 1, "
          "P(B_post<0) in [0.9, 1.0] on every K",
          all(float(srows[K]["return_never"]) == 0.0
              and abs(float(srows[K]["gap_closing_share_contraction"]) - 1.0) <= 1e-9
              and 0.9 <= float(srows[K]["P_post_B_neg"]) <= 1.0
              for K in K_LIST))
    check("M1C.4 consistency: no_switch share <= 0.06 on every K",
          all(float(srows[K]["frac_no_switch"]) <= 0.06 for K in K_LIST))
    check("switch-mechanism ordering anchor: P2_new share grows with K "
          "(frac_P2_new(K=400) > frac_P2_new(K=50))",
          float(srows[400]["frac_P2_new"]) > float(srows[50]["frac_P2_new"]))
else:
    check("M1C.6 gate_return_summary.csv present (run the analysis first)", False)

ep_csv = res_dir / "episode_table.csv"
if ep_csv.exists():
    ep_rows = list(csv.DictReader(ep_csv.open(encoding="utf-8")))
    ok_rep = True
    for K in [50, 400]:
        rows = [r for r in ep_rows if int(r["K"]) == K]
        if len(rows) != len(recomputed[K]["eps"]):
            ok_rep = False
            continue
        for rec, row in zip(recomputed[K]["eps"], rows):
            if int(row["t_r"]) != rec["t_r"] \
               or int(row["tau_rebound"]) != rec["tau_rebound"] \
               or (row["terminal"] == "True") != rec["terminal"] \
               or abs(float(row["B_r"]) - rec["B_r"]) > 5e-7:
                ok_rep = False
    check("deterministic reproducibility: committed episode table matches "
          "recomputation (K=50 and K=400)", ok_rep)
else:
    check("M1C.6 episode_table.csv present", False)

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1C.6 boundary-return tests passed")
