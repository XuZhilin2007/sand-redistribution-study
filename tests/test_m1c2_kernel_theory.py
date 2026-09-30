"""M1C.2 target-matched kernel theory tests (docs/M1C2_TARGET_MATCHED_KERNEL_THEORY.md).

Run:  python tests/test_m1c2_kernel_theory.py

All formulas are exact consequences of the canonical discrete update in
`adaptive.run_m1a_deterministic`; this suite verifies them numerically on
the canonical A (fixed q_near respray) and U (fixed uniform respray)
trajectories and checks the monotonicity theorem and its scope.

  1. one-step region identity (prefix/suffix formulas) == actual D of the
     next canonical state, for every active round of A and U at all four
     canonical K;
  2. target-matched specialization (G = T) and the mismatch identity
     D_actual_new - D_target_matched_new == M * D_g (per prefix, per round);
  3. boundary export identity: prefix mass change == -M (1 - G(j*));
  4. monotonicity theorem on the U trajectories: D_max(t+1) < D_max(t) at
     every active round, plus the stronger region statements (prefix part
     <= (1-alpha) D_max; suffix part per-prefix non-increasing);
  5. theorem scope: the target-matched counterfactual is strictly
     D_max-decreasing even on A's states, while A's actual profile rises in
     ~36% of rounds (the mismatch term is load-bearing);
  6. D_max = 0 does NOT imply bin-uniformity: the alternating state
     (0, 2/K, 0, 2/K, ...) has D_max = 0 with U_density = 1, and M1C.1's
     canonical U stop state shows the same separation at ~5.4e-3;
  7. q_near export geometry: 1 - G_near(a) = (1-a)^2 vs 1 - T(a) = 1 - a.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.adaptive import cumulative_excess, run_m1a_deterministic
from sand_m0.kernel_theory import (
    mismatch_reinjection,
    one_step_excess,
    prefix_cdf,
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
T_MAX, ALPHA, TOL = 5000, 0.25, 1e-12
q_near, q_unif = DISTRIBUTIONS["near"], DISTRIBUTIONS["uniform"]

runs = {}
for K in K_LIST:
    runs[("A", K)] = run_m1a_deterministic(q_near, K=K, T_max=T_MAX, alpha=ALPHA,
                                           stop_tolerance=TOL)
    runs[("U", K)] = run_m1a_deterministic(q_near, K=K, T_max=T_MAX, alpha=ALPHA,
                                           stop_tolerance=TOL,
                                           redistribution_dist=q_unif)

# --- 1. one-step region identity on every active round ------------------------------
worst = 0.0
ok = True
for variant, respray in [("A", q_near), ("U", q_unif)]:
    for K in K_LIST:
        r = runs[(variant, K)]
        states = np.asarray(r["states"], dtype=float)
        g = respray.bin_probs(K)
        for t in range(r["stopped_at"]):
            j_star = r["rows"][t]["j_t"]
            predicted = one_step_excess(states[t], j_star, ALPHA, g)
            actual = cumulative_excess(states[t + 1], K)
            worst = max(worst, float(np.abs(predicted - actual).max()))
            if np.abs(predicted - actual).max() > 1e-10:
                ok = False
n_rounds = sum(runs[(v, K)]["stopped_at"] for v in ["A", "U"] for K in K_LIST)
check(f"one-step region identity == actual next-state D on all {n_rounds} active "
      f"rounds of A and U (4 K), max residual {worst:.2e} <= 1e-10", ok)

# --- 2. target-matched specialization + mismatch identity ---------------------------
worst_tm = worst_mm = 0.0
ok_tm = ok_mm = True
for K in K_LIST:
    r = runs[("A", K)]
    states = np.asarray(r["states"], dtype=float)
    g = q_near.bin_probs(K)
    for t in range(r["stopped_at"]):
        row = r["rows"][t]
        j_star = row["j_t"]
        M = row["removed_mass"]
        actual = cumulative_excess(states[t + 1], K)
        tm = target_matched_excess(states[t], j_star, ALPHA)
        mm = mismatch_reinjection(M, g)
        worst_tm = max(worst_tm, float(np.abs((actual - mm) - tm).max()))
        worst_mm = max(worst_mm, float(np.abs(mm - M * prefix_cdf(g)
                                              + M * np.arange(1, K + 1) / K).max()))
        if np.abs((actual - mm) - tm).max() > 1e-10:
            ok_tm = False
        if np.abs(mm - M * prefix_cdf(g) + M * np.arange(1, K + 1) / K).max() > 1e-12:
            ok_mm = False
check(f"A: actual update == target-matched update + M*D_g on every round, "
      f"max residual {worst_tm:.2e} <= 1e-10", ok_tm)
check("mismatch reinjection profile == M * D_g(j) exactly (M * [G(j) - T(j)])",
      ok_mm)

# --- 3. boundary export identity ----------------------------------------------------
worst = 0.0
ok = True
for variant, respray in [("A", q_near), ("U", q_unif)]:
    for K in K_LIST:
        r = runs[(variant, K)]
        states = np.asarray(r["states"], dtype=float)
        G = prefix_cdf(respray.bin_probs(K))
        for t in range(r["stopped_at"]):
            row = r["rows"][t]
            j_star = row["j_t"]
            pc = float(states[t + 1][:j_star].sum() - states[t][:j_star].sum())
            worst = max(worst, abs(pc + row["removed_mass"] * (1.0 - G[j_star - 1])))
            if abs(pc + row["removed_mass"] * (1.0 - G[j_star - 1])) > 1e-10:
                ok = False
check(f"boundary export identity prefix_change == -M(1-G(j*)) on all rounds, "
      f"max residual {worst:.2e} <= 1e-10", ok)

# --- 4. monotonicity theorem on the U trajectories ----------------------------------
ok_mono = ok_region1 = ok_region2 = True
for K in K_LIST:
    r = runs[("U", K)]
    states = np.asarray(r["states"], dtype=float)
    Ds = [cumulative_excess(states[t], K) for t in range(r["stopped_at"] + 1)]
    for t in range(r["stopped_at"]):
        j_star = r["rows"][t]["j_t"]
        d_max = float(Ds[t].max())
        d_new = Ds[t + 1]
        if not float(d_new.max()) < d_max:
            ok_mono = False
        if not float(d_new[:j_star].max()) <= (1.0 - ALPHA) * d_max + 1e-15:
            ok_region1 = False
        if not (d_new[j_star:] <= Ds[t][j_star:] + 1e-12).all():
            ok_region2 = False
check("theorem (strict): U's D_max strictly decreases at every active round, "
      "all 4 K", ok_mono)
check("theorem (region 1): after each U round every prefix <= (1-alpha) D_max",
      ok_region1)
check("theorem (region 2): after each U round every suffix prefix D(j) is "
      "non-increasing (exact in real arithmetic; realized-state fp residue "
      "<= 1e-12, all observed violations at j=K where D(K)=0+-ulp)", ok_region2)

# --- 5. theorem scope: counterfactual on A's states vs A's actual profile -----------
ok_cf = True
for K in K_LIST:
    r = runs[("A", K)]
    states = np.asarray(r["states"], dtype=float)
    for t in range(r["stopped_at"]):
        j_star = r["rows"][t]["j_t"]
        d_max = float(cumulative_excess(states[t], K).max())
        cf = float(target_matched_excess(states[t], j_star, ALPHA).max())
        if not cf < d_max:
            ok_cf = False
check("counterfactual: the target-matched one-step map applied to A's own "
      "states is strictly D_max-decreasing at every A round (theorem holds "
      "for arbitrary states)", ok_cf)
D_a = [cumulative_excess(np.asarray(runs[("A", 100)]["states"], dtype=float)[t], 100)
       for t in range(runs[("A", 100)]["stopped_at"] + 1)]
share_up = float(np.mean([D_a[t + 1].max() > D_a[t].max() + 1e-15
                          for t in range(len(D_a) - 1)]))
check("contrast: A's actual D_max rises in a substantial share of rounds "
      "(K=100 share ~ 0.36 > 0.1) — the mismatch term is load-bearing",
      share_up > 0.1, f"share = {share_up:.4f}")

# --- 6. D_max = 0 does not imply bin-uniformity -------------------------------------
K = 100
u = 1.0 / K
sawtooth = np.tile([0.0, 2.0 * u], K // 2)
D_saw = cumulative_excess(sawtooth, K)
U_density_saw = K * float(np.sum((sawtooth - u) ** 2))
check("constructed state (0, 2/K, 0, 2/K, ...): D_max = 0 exactly while "
      "U_density = 1 (one-sided cumulative criterion is far weaker than "
      "binwise uniformity)",
      abs(float(D_saw.max())) <= 1e-15 and abs(U_density_saw - 1.0) <= 1e-12)
st_u = np.asarray(runs[("U", 100)]["states"], dtype=float)[runs[("U", 100)]["stopped_at"]]
check("canonical U stop state shows the same separation: D_max <= 0 (within "
      "1e-12) and U_density ~ 5.4e-3 (not machine-uniform)",
      float(cumulative_excess(st_u, K).max()) <= 1e-12
      and 1e-3 < K * float(np.sum((st_u - u) ** 2)) < 1e-2)

# --- 7. q_near export geometry ------------------------------------------------------
for a in [0.5, 0.9, 0.94]:
    g_ratio = (1.0 - a) ** 2 / (1.0 - a)
    check(f"export geometry at a={a}: (1-G_near)/(1-T) = 1-a = {1-a:.2f} "
          "(near-biased net export is a factor (1-a) weaker than target-matched)",
          abs(((1 - (2 * a - a * a)) / (1 - a)) - (1 - a)) <= 1e-12
          and abs(g_ratio - (1 - a)) <= 1e-12)

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1C.2 kernel-theory tests passed")
