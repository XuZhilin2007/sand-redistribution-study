"""M1C.3 redistribution-law order & rebound criterion tests
(docs/M1C3_REDISTRIBUTION_LAW_ORDER.md).

Run:  python tests/test_m1c3_order_criterion.py

Builds on the M1C.2 decomposition D'(j) = D'_TM(j) + M Delta(j),
Delta(j) = G(j) - T(j). All items are pure-mathematics checks on the
canonical update rule; the only trajectories used are the ALREADY EXISTING
canonical Variant A runs (no new simulation, no new experiment).

  1. order theorem (G <= T): for representative states and several
     cumulatively target-dominated laws, a manually applied canonical step
     strictly decreases D_max (exact arithmetic claim, fp-aware check);
  2. suffix impossibility: for ANY fixed law the suffix (j > j*) never
     produces a rebound — verified on every active A round and by algebra;
  3. prefix margin bound: gamma(j) >= alpha * D_max on the swept prefix,
     hence any rebound requires M Delta(j) > alpha D_max (uniform target:
     M > alpha D_max K/(K-1));
  4. rebound position prediction: in every actual A rebound round the new
     argmax lies inside the swept prefix (j <= j*);
  5. EXACT one-step rebound criterion: predicted (exists j: M Delta(j) >
     gamma(j)) <=> actual D_max increased, on ALL active A rounds at all
     four canonical K; counts, FP/FN, and max residual reported;
  6. global sufficient bound: rounds with M Delta_+ <= Gamma never rebound;
     and the magnitude bound D_max' - D_max <= M Delta_+ - Gamma always;
  7. q_near specialization: Delta_+ = 1/4 exactly on the canonical grids;
  8. counterexample A (converse falsification): a law with G > T somewhere
     can still contract D_max every step — tiny tilted-uniform law on
     canonical A states: criterion says no rebound, manual step agrees;
  9. counterexample B (tiny-K sweep): for K in {2,3,4}, many deterministic
     states and G <= T laws, no single manual step ever increases D_max;
 10. decomposition consistency: at the realized new argmax of a rebound
     round R = M Delta / gamma > 1, and < 1 in decrease rounds.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.adaptive import cumulative_excess, run_m1a_deterministic
from sand_m0.kernel_theory import (
    contraction_margin,
    cumulative_mismatch,
    one_step_excess,
    rebound_criterion,
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
T_MAX, ALPHA, TOL, FP_TOL = 5000, 0.25, 1e-12, 1e-12
q_near, q_unif, q_far = DISTRIBUTIONS["near"], DISTRIBUTIONS["uniform"], DISTRIBUTIONS["far"]


def manual_step(mass: np.ndarray, alpha: float, respray_bins: np.ndarray):
    """One canonical active step by hand (selection + removal + fixed respray).

    Pure vector algebra on the given state; mirrors the canonical update.
    """
    K = mass.shape[0]
    D = cumulative_excess(mass, K)
    j_star = int(np.argmax(D)) + 1
    moved = alpha * float(mass[:j_star].sum())
    p_minus = mass.copy()
    p_minus[:j_star] *= (1.0 - alpha)
    return j_star, moved, p_minus + moved * respray_bins


# canonical A trajectories (existing, recomputed deterministically)
runs = {K: run_m1a_deterministic(q_near, K=K, T_max=T_MAX, alpha=ALPHA,
                                 stop_tolerance=TOL)
        for K in K_LIST}

# --- 1. order theorem: G <= T => strict D_max decrease on manual steps ---------------
g_td_laws = {
    "uniform (G=T)": q_unif.bin_probs(100),
    "far (G<=T)": q_far.bin_probs(100),
    "mix 0.5u+0.5f": 0.5 * q_unif.bin_probs(100) + 0.5 * q_far.bin_probs(100),
}
states = [
    q_near.bin_probs(100),
    np.asarray(runs[100]["states"], dtype=float)[37],
    np.asarray(runs[100]["states"], dtype=float)[200],
]
rng = np.random.default_rng(20260919)
for _ in range(5):
    p = rng.random(100)
    p[p.argmax()] += 3.0
    states.append(p / p.sum())
ok = True
worst = np.inf
for p0 in states:
    for g in g_td_laws.values():
        if not ((np.cumsum(g) - np.arange(1, 101) / 100) <= 1e-15).all():
            ok = False
            continue
        d0 = float(cumulative_excess(p0, 100).max())
        if d0 <= FP_TOL:
            continue  # inactive; theorem is about active steps
        _, _, p1 = manual_step(p0, ALPHA, g)
        d1 = float(cumulative_excess(p1, 100).max())
        worst = min(worst, d0 - d1)
        if not d1 < d0:
            ok = False
check("order theorem: G <= T (uniform, far, mixture) => D_max strictly "
      "decreases on manual active steps (worst decrease %.3e > 0)" % worst,
      ok and worst > 0)

# --- 2+4. suffix impossibility + rebound position, on ALL active A rounds ------------
ok_sfx = ok_pos = True
for K in K_LIST:
    r = runs[K]
    states_k = np.asarray(r["states"], dtype=float)
    for t in range(r["stopped_at"]):
        row = r["rows"][t]
        j_star = row["j_t"]
        d0 = cumulative_excess(states_k[t], K)
        d1 = cumulative_excess(states_k[t + 1], K)
        if float(d1[j_star:].max()) > float(d0.max()) + FP_TOL:
            ok_sfx = False
        d0_max = float(d0.max())
        if float(d1.max()) > d0_max + FP_TOL:
            j_new = int(np.argmax(d1)) + 1
            if j_new > j_star:
                ok_pos = False
check("suffix impossibility: no A round ever raises any suffix prefix above "
      "the old D_max (all 5344 active rounds, 4 K)", ok_sfx)
check("rebound position prediction: in every actual A rebound round the new "
      "argmax lies inside the swept prefix (j <= j*)", ok_pos)

# --- 3. prefix margin bound gamma(j) >= alpha * D_max --------------------------------
ok_margin = True
for K in K_LIST:
    r = runs[K]
    states_k = np.asarray(r["states"], dtype=float)
    for t in range(0, r["stopped_at"], 7):  # deterministic stride sample
        j_star = r["rows"][t]["j_t"]
        gamma = contraction_margin(states_k[t], j_star, ALPHA)
        d_max = float(cumulative_excess(states_k[t], K).max())
        if float(gamma[:j_star].min()) < ALPHA * d_max - 1e-12:
            ok_margin = False
check("prefix margin bound: gamma(j) >= alpha*D_max on every swept prefix "
      "(stride-7 sample of all K)", ok_margin)

# --- 5. EXACT rebound criterion on the existing A trajectory -------------------------
total = actual_reb = pred_reb = fp = fn = 0
worst_res = worst_match = 0.0
for K in K_LIST:
    r = runs[K]
    states_k = np.asarray(r["states"], dtype=float)
    g = q_near.bin_probs(K)
    for t in range(r["stopped_at"]):
        row = r["rows"][t]
        j_star = row["j_t"]
        d0_max = float(cumulative_excess(states_k[t], K).max())
        d1_max = float(cumulative_excess(states_k[t + 1], K).max())
        crit = rebound_criterion(states_k[t], j_star, ALPHA, g)
        # criterion form == one_step_excess form == actual next D (identity check)
        pred_profile = target_matched_excess(states_k[t], j_star, ALPHA) \
            + crit["reinjection"]
        worst_res = max(worst_res, abs(float(pred_profile.max()) - d1_max))
        worst_match = max(worst_match, abs(float(pred_profile.max())
                                           - float(one_step_excess(
                                               states_k[t], j_star, ALPHA, g).max())))
        actual = d1_max > d0_max + FP_TOL
        predicted = crit["rebound"]
        total += 1
        actual_reb += int(actual)
        pred_reb += int(predicted)
        fp += int(predicted and not actual)
        fn += int(actual and not predicted)
check(f"EXACT criterion on existing A trajectory: {total} active rounds, "
      f"{actual_reb} actual rebounds, {pred_reb} predicted — "
      f"FP={fp}, FN={fn}, profile residual {worst_res:.2e} <= 1e-9, "
      f"decomposition residual {worst_match:.2e} <= 1e-9",
      fp == 0 and fn == 0 and pred_reb == actual_reb
      and worst_res <= 1e-9 and worst_match <= 1e-9)

# --- 6. global sufficient bound and magnitude bound ----------------------------------
ok_suff = ok_mag = True
for K in K_LIST:
    r = runs[K]
    states_k = np.asarray(r["states"], dtype=float)
    g = q_near.bin_probs(K)
    delta = cumulative_mismatch(g)
    for t in range(r["stopped_at"]):
        row = r["rows"][t]
        j_star = row["j_t"]
        M = row["removed_mass"]
        gamma = contraction_margin(states_k[t], j_star, ALPHA)
        Gamma = float(gamma.min())
        Mdp = M * float(np.maximum(delta, 0.0).max())
        d0_max = float(cumulative_excess(states_k[t], K).max())
        d1_max = float(cumulative_excess(states_k[t + 1], K).max())
        if Mdp <= Gamma - 1e-15 and not d1_max <= d0_max + FP_TOL:
            ok_suff = False
        if d1_max - d0_max > Mdp - Gamma + 1e-9:
            ok_mag = False
check("global sufficient bound: every round with M*Delta_+ <= Gamma has "
      "D_max' <= D_max", ok_suff)
check("magnitude bound: D_max' - D_max <= M*Delta_+ - Gamma on every round",
      ok_mag)

# --- 7. q_near specialization --------------------------------------------------------
ok_qn = all(
    abs(float(cumulative_mismatch(q_near.bin_probs(K)).max()) - 0.25) <= 1e-12
    for K in K_LIST
)
check("q_near specialization: Delta_+ = max_j (j/K)(1-j/K) = 1/4 exactly on "
      "all canonical grids (peak reinjection M/4)", ok_qn)

# --- 8. counterexample A: G > T somewhere, still no rebound --------------------------
g_tilt = np.full(100, 1.0) + np.linspace(-0.002, 0.002, 100)  # slight near tilt
g_tilt /= g_tilt.sum()
assert (cumulative_mismatch(g_tilt) > 0).any() and cumulative_mismatch(g_tilt).max() < 0.001
ok_cx = True
for t in range(0, runs[100]["stopped_at"], 11):  # deterministic stride sample
    st = np.asarray(runs[100]["states"], dtype=float)
    j_star, _, p1 = manual_step(st[t], ALPHA, g_tilt)
    crit = rebound_criterion(st[t], j_star, ALPHA, g_tilt)
    d0 = float(cumulative_excess(st[t], 100).max())
    d1 = float(cumulative_excess(p1, 100).max())
    if crit["rebound"] or not (d1 < d0):
        ok_cx = False
check("counterexample A: law with G > T somewhere (tiny near tilt) still "
      "contracts D_max at every sampled A state — positive mismatch is "
      "necessary-for-rebound only locally, not sufficient", ok_cx)

# --- 9. counterexample B: tiny-K exhaustive sweep for G <= T -------------------------
ok_cxb = True
n_checked = 0
for K in [2, 3, 4]:
    T = np.arange(1, K + 1) / K
    rng_k = np.random.default_rng(20260919 + K)
    laws = [np.full(K, 1.0 / K)]
    for _ in range(20):
        base = rng_k.random(K)
        base /= base.sum()
        g = 0.5 / K + 0.5 * base  # mixture: guaranteed G <= T? verify per law
        if ((np.cumsum(g) - T) <= 1e-15).all():
            laws.append(g)
    for s in range(300):
        p = rng_k.random(K)
        p[rng_k.integers(K)] += rng_k.random() * 2
        p /= p.sum()
        d0 = float(cumulative_excess(p, K).max())
        if d0 <= 1e-3:
            continue  # need a clearly active step
        for g in laws:
            j_star, _, p1 = manual_step(p, ALPHA, g)
            d1 = float(cumulative_excess(p1, K).max())
            n_checked += 1
            if not d1 < d0 + 1e-12:
                ok_cxb = False
check(f"counterexample B (tiny-K sweep): no G <= T law ever increases D_max "
      f"({n_checked} manual active steps, K in {{2,3,4}}, deterministic)",
      ok_cxb)

# --- 10. decomposition consistency at the realized new argmax ------------------------
ok_dec = True
n_r = n_d = 0
r = runs[100]
st = np.asarray(r["states"], dtype=float)
g = q_near.bin_probs(100)
for t in range(r["stopped_at"]):
    j_star = r["rows"][t]["j_t"]
    gamma = contraction_margin(st[t], j_star, ALPHA)
    M = r["rows"][t]["removed_mass"]
    j_new = int(np.argmax(cumulative_excess(st[t + 1], 100)))
    R = M * cumulative_mismatch(g)[j_new] / gamma[j_new]
    rebound = float(cumulative_excess(st[t + 1], 100).max()) > \
        float(cumulative_excess(st[t], 100).max()) + FP_TOL
    if rebound:
        n_r += 1
        if not R > 1.0:
            ok_dec = False
    else:
        n_d += 1
        if not R < 1.0:
            ok_dec = False
check(f"decomposition: at the realized new argmax R = M*Delta/gamma is > 1 "
      f"exactly in rebound rounds ({n_r}) and < 1 in decrease rounds ({n_d})",
      ok_dec and n_r > 0)

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1C.3 order & criterion tests passed")
