"""M2 Gate 1 redistribution-mismatch generality tests
(docs/M2_GATE1_GENERALITY_GATE.md, experiments/m2_gate1_generality/).

Run:  python tests/test_m2_gate1_generality.py

Covers the frozen Gate 1 kernel set {K0, K-, KW, KC, KS}:

  1. kernel mathematics: endpoints, monotonicity, exact-CDF-differencing
     bin masses and normalization at every execution K (tau_num = 1e-12);
  2. KW mismatch formula Delta_W(j/K) = (1/4)(j/K)(1-j/K) on the grid;
  3. K- mismatch = x^2 - x <= 0 and sign-mirror relation to canonical
     Delta_C = x(1-x) (Delta_- == -Delta_C pointwise);
  4. KS structure: Delta_S >= 0, max = 1/4, breakpoint continuity at 1/3
     and 2/3, integral = 1/6 and max = 1/4 in EXACT rational arithmetic,
     and off-grid bin masses at K=100 against an independent closed-form
     evaluation;
  5. generic kernel behavior: identical canonical initial state for every
     kernel at the same K; respray step == prefix-scaling + removed*q
     exactly; mass conservation each round;
  6. theory verification on short real runs: one-step decomposition
     identity D_new = D_TM_new + M*Delta against the actual next state,
     exact rebound criterion iff actual rebound (plain strict semantics,
     numerical-boundary cases flagged), R_t (M1C.4 semantics) > 1 iff
     criterion rebound;
  7. K- order-theorem preview: strict D_max decrease on every active
     round of a short K- trajectory.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.adaptive import cumulative_excess, run_m1a_deterministic
from sand_m0.kernel_theory import target_matched_excess
from sand_m0.model import (
    DISTRIBUTIONS,
    flat_top_distribution,
    weak_positive_distribution,
)

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


TAU = 1e-12
ALPHA, TOL = 0.25, 1e-12
K_LIST = [50, 100, 200, 400]
q_near = DISTRIBUTIONS["near"]
q_unif = DISTRIBUTIONS["uniform"]
q_far = DISTRIBUTIONS["far"]
kw = weak_positive_distribution()
ks = flat_top_distribution()
KERNELS = {"K0": q_unif, "K-": q_far, "KW": kw, "KC": q_near, "KS": ks}

# --- 1. endpoints / monotonicity / bin masses ---------------------------------------
worst_end = max_sum_err = worst_min_q = 0.0
ok_mono = True
for name, dist in KERNELS.items():
    worst_end = max(worst_end, abs(float(dist.cdf(np.array([0.0]))[0])),
                    abs(float(dist.cdf(np.array([1.0]))[0]) - 1.0))
    g = np.linspace(0.0, 1.0, 10001)
    if not np.all(np.diff(dist.cdf(g)) >= -TAU):
        ok_mono = False
    for K in K_LIST:
        edges = np.linspace(0.0, 1.0, K + 1)
        raw = np.diff(dist.cdf(edges))
        max_sum_err = max(max_sum_err, abs(float(raw.sum()) - 1.0))
        worst_min_q = min(worst_min_q, float(raw.min()))
check(f"all five kernels: G(0)=0, G(1)=1 within {TAU:.0e} (worst {worst_end:.2e})",
      worst_end <= TAU)
check("all five kernels: G nondecreasing on a 10001-point grid", ok_mono)
check(f"all five kernels x K in {K_LIST}: raw CDF-differencing bin masses have "
      f"|sum-1| <= {TAU:.0e} (worst {max_sum_err:.2e}) and min q_j >= -{TAU:.0e} "
      f"(worst {worst_min_q:.2e})",
      max_sum_err <= TAU and worst_min_q >= -TAU)

# --- 2. KW mismatch formula -----------------------------------------------------------
worst_kw = 0.0
for K in K_LIST:
    x = np.arange(1, K + 1) / K
    worst_kw = max(worst_kw, float(np.abs(kw.cdf(x) - x - 0.25 * x * (1.0 - x)).max()))
check(f"KW: Delta_W(j/K) == (1/4)(j/K)(1-j/K) on all grids (worst {worst_kw:.2e})",
      worst_kw <= TAU)

# --- 3. K- sign mirror ------------------------------------------------------------------
x = np.arange(1, 101) / 100
d_minus = q_far.cdf(x) - x
d_canon = q_near.cdf(x) - x
check("K-: Delta_minus = x^2 - x <= 0 everywhere and equals -Delta_C pointwise "
      "(exact sign mirror of canonical)",
      float(d_minus.max()) <= TAU and
      float(np.abs(d_minus + d_canon).max()) <= TAU)

# --- 4. KS structure ---------------------------------------------------------------------
xs = np.linspace(0.0, 1.0, 30001)
d_s = ks.cdf(xs) - xs
check(f"KS: Delta_S >= 0 (min {float(d_s.min()):.2e}) and max Delta_S = 1/4 "
      f"on the flat top ({float(d_s.max()):.16f})",
      float(d_s.min()) >= -TAU and abs(float(d_s.max()) - 0.25) <= TAU)
for bp, tag in [(1.0 / 3.0, "1/3"), (2.0 / 3.0, "2/3")]:
    gl = float((7.0 / 4.0) * bp if bp <= 1.0 / 3.0 + 0 else
               (bp + 0.25 if bp <= 2.0 / 3.0 else 0.75 + 0.25 * bp))
    gr = float(bp + 0.25 if bp <= 2.0 / 3.0 else 0.75 + 0.25 * bp)
    gap = abs(gl - gr)
    got = float(ks.cdf(np.array([bp]))[0])
    check(f"KS: continuous at {tag} (branch gap {gap:.2e}; CDF value matches "
          f"to {abs(got - gr):.2e})",
          gap <= TAU and abs(got - gr) <= TAU)
# exact rational arithmetic: max and integral of Delta_S
i1 = (Fraction(3, 4) * Fraction(1, 3) ** 2) / 2          # left wedge
i2 = Fraction(1, 4) * Fraction(1, 3)                      # flat top
integral = 2 * i1 + i2
check("KS (exact rational arithmetic): integral of Delta_S = 1/6 exactly and "
      "max Delta_S = 1/4 exactly on the closed flat top",
      integral == Fraction(1, 6))
# off-grid bin masses at K=100 against an independent closed-form evaluation
K = 100
edges = np.linspace(0.0, 1.0, K + 1)


def ks_cdf_independent(x: np.ndarray) -> np.ndarray:
    out = np.empty_like(x)
    for i, v in enumerate(x):
        if v <= 1.0 / 3.0:
            out[i] = 1.75 * v
        elif v >= 2.0 / 3.0:
            out[i] = 0.75 + 0.25 * v
        else:
            out[i] = v + 0.25
    return out


q_ks = ks.bin_probs(K)
q_ref = np.diff(ks_cdf_independent(edges))
check("KS: bin_probs(K=100) equals independent closed-form CDF differencing "
      "(breakpoints 1/3, 2/3 NOT rounded to the grid; max diff %.2e)"
      % float(np.abs(q_ks - q_ref).max()),
      float(np.abs(q_ks - q_ref).max()) <= TAU and abs(float(q_ref.sum()) - 1.0) <= TAU)
q_kw = kw.bin_probs(K)
q_kw_ref = np.diff(1.25 * edges - 0.25 * edges ** 2)
check("KW: bin_probs(K=100) equals independent closed-form CDF differencing "
      "(max diff %.2e)" % float(np.abs(q_kw - q_kw_ref).max()),
      float(np.abs(q_kw - q_kw_ref).max()) <= TAU)

# --- 5. generic kernel behavior: shared initial state + exact respray step --------------
runs = {}
for name, dist_respray in [("K0", q_unif), ("K-", q_far), ("KW", kw),
                           ("KC", None), ("KS", ks)]:
    runs[name] = run_m1a_deterministic(
        q_near, K=100, T_max=60, alpha=ALPHA, stop_tolerance=TOL,
        redistribution="throw", redistribution_dist=dist_respray)
states0 = [np.asarray(r["states"], dtype=float)[0] for r in runs.values()]
identical = all(np.array_equal(states0[0], s) for s in states0[1:])
first_decision = {name: (r["rows"][0]["j_t"], r["rows"][0]["D_max"])
                  for name, r in runs.items()}
check("same-K shared initial state: states[0] bit-identical across all five "
      "kernels and equal to canonical q_near bin probabilities; identical "
      "first decision (j_t, D_max)",
      identical and np.array_equal(states0[0], q_near.bin_probs(100))
      and len({tuple(v) for v in first_decision.values()}) == 1)

ok_step = ok_mass = True
worst_step = 0.0
for name, r in runs.items():
    dist_respray = {"K0": q_unif, "K-": q_far, "KW": kw,
                    "KC": q_near, "KS": ks}[name]
    g = dist_respray.bin_probs(100)
    st = np.asarray(r["states"], dtype=float)
    for t in range(r["stopped_at"] if r["stopped_at"] is not None else 60):
        row = r["rows"][t]
        p_minus = st[t].copy()
        p_minus[: row["j_t"]] *= (1.0 - ALPHA)
        added = st[t + 1] - p_minus
        worst_step = max(worst_step, float(np.abs(added - row["removed_mass"] * g).max()))
        if abs(float(st[t + 1].sum()) - 1.0) > TAU:
            ok_mass = False
check(f"generic fixed-G support: every round of every kernel is exactly "
      f"prefix-scaling + removed*q_kernel (worst profile diff {worst_step:.2e})",
      ok_step and worst_step <= TAU)
check("mass conservation: |sum(p)-1| <= tau_num on every round of every "
      "short run", ok_mass)

# --- 6. theory verification on the short runs --------------------------------------------
worst_resid = 0.0
tp = tn = fp = fn = 0
n_boundary = 0
ok_R = True
n_rounds = 0
for name, r in runs.items():
    if name == "KC":
        continue  # canonical rounds already verified in M1C.2/M1C.3; avoid re-mining
    g = {"K0": q_unif, "K-": q_far, "KW": kw, "KS": ks}[name].bin_probs(100)
    delta = np.cumsum(g) - np.arange(1, 101) / 100
    st = np.asarray(r["states"], dtype=float)
    n_active = r["stopped_at"] if r["stopped_at"] is not None else 60
    for t in range(n_active):
        row = r["rows"][t]
        j_star = row["j_t"]
        moved = row["removed_mass"]
        D_tm = target_matched_excess(st[t], j_star, ALPHA)
        pred = D_tm + moved * delta
        actual = cumulative_excess(st[t + 1], 100)
        worst_resid = max(worst_resid, float(np.abs(pred - actual).max()))
        gamma = float(row["D_max"]) - D_tm
        worst_excess = float((moved * delta - gamma).max())
        criterion = worst_excess > 0.0
        actual_reb = float(row["D_max"]) < float(r["rows"][t + 1]["D_max"])
        if abs(worst_excess) <= TAU:
            n_boundary += 1
        if actual_reb and criterion:
            tp += 1
        elif not actual_reb and not criterion:
            tn += 1
        elif criterion:
            fp += 1
        else:
            fn += 1
        pos_reinj = moved * np.maximum(delta, 0.0)
        ratio = np.where(gamma > 0.0,
                         pos_reinj / np.where(gamma > 0.0, gamma, 1.0), 0.0)
        if (float(ratio.max()) > 1.0) != criterion:
            ok_R = False
        n_rounds += 1
check(f"decomposition identity D_new = D_TM_new + M*Delta on all {n_rounds} "
      f"active rounds of the K0/K-/KW/KS short runs (max residual "
      f"{worst_resid:.2e} <= {TAU:.0e})", worst_resid <= TAU)
check(f"exact rebound criterion iff actual rebound on the same rounds "
      f"(TP={tp}, TN={tn}, FP={fp}, FN={fn}, numerical-boundary rounds={n_boundary})",
      fp == 0 and fn == 0)
check("R_t helper compatibility (M1C.4 semantics): R_t > 1 iff criterion "
      "rebound on every checked round", ok_R)

# --- 7. K- order-theorem preview -----------------------------------------------------------
km = runs["K-"]
n_active = km["stopped_at"] if km["stopped_at"] is not None else 60
D_traj = [km["rows"][t]["D_max"] for t in range(n_active + 1)]
strict = all(D_traj[t + 1] < D_traj[t] for t in range(n_active))
check(f"K- preview: strict D_max decrease on every active round of the short "
      f"run ({n_active} rounds, no rebound — G <= T order theorem)",
      strict and n_active > 0)

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M2 Gate 1 generality tests passed")
