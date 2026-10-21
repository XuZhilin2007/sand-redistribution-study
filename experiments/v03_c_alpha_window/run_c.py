"""v0.3 Theorem C runner: alpha-universal finite-window classification.

Exact verification on an (alpha, K) grid that:
  C1  dwell bounds actually hold on the exact window itineraries
      (entry t* <= T_alpha; A-run lengths <= max(N_A^alpha, U1^alpha);
      B-run lengths <= N_B^alpha);
  C2  the two-test classifier (any ell-selection after tau_0 within the
      window <=> CYCLE) is exactly correct with window
      H(alpha,K) = floor(81K/(80 alpha)) + max(N_B, U1) + N_A + 2;
  C3  the frozen priority classifier (LOCK iff trailing m-run >= K) is
      guaranteed correct when N_A^alpha < K (checked at the longer window
      H'), and the registered counterexample K=9, alpha=1/8
      (N_A = 18 >= 9) misclassifies LOCK at its registered window H=90
      while the two-test classifier is correct;
  C4  LOCK classes never select ell after t* + N_B^alpha (within window).

For LOCK classes (L >= h) the A-run bounds N_A/U1 are defined by the
uniform analytic bound min{n: (1-alpha)^n <= 1/81} (the crossing-based
loops diverge there); they only keep the window formula finite.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import sys
import time
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "s4_w1_s1_xo_precision"))

import replay  # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

CFG = {
    "alphas": [F(1, 4), F(1, 8), F(1, 100), F(1, 2), F(3, 4), F(2, 3)],
    "Ks": [9, 10, 11, 15, 16, 17, 21, 22, 23, 33, 34, 35, 51, 52, 53,
           12, 13, 14, 24, 25, 26],
    "Ks_cap_small_alpha": [9, 10, 11, 15, 16, 17, 21, 22, 23,
                           12, 13, 14],
    "counterexample": {"K": 9, "alpha": F(1, 8), "H": 90},
}


def params(K):
    m = (5 * K + 5) // 6  # == ceil(5K/6), pure integer arithmetic
    h = F(1, K)
    b = F(m, K)
    c = 2 * b - b * b
    q = F(6 - K % 6, 3 * K)
    L = c * q
    return m, h, b, c, q, L


def bounds(K, alpha):
    m, h, b, c, q, L = params(K)
    one_minus = 1 - alpha
    if L >= h:
        # LOCK class: post-landing/first A-runs never end (x stays above
        # h); the A-run bounds are vacuous. Use the uniform analytic bound
        # (1-alpha)^n <= 1/81 (lambda <= 80/81 across all K) so that the
        # window formula stays finite.
        NA = 1
        while one_minus ** NA > F(1, 81):
            NA += 1
        U1 = NA
    else:
        # post-landing A-run: (1-alpha)^n (h + alpha L - L) <= h - L
        NA = 1
        while one_minus ** NA * (h - (1 - alpha) * L) > h - L:
            NA += 1
        # first A-run from x0 <= c: (1-alpha)^n (c - L) <= h - L
        U1 = 1
        while one_minus ** U1 * (c - L) > h - L:
            U1 += 1
    # initial B-run: (1 - alpha q)^n (c) < c - h  (worst x0 = 0)
    f = 1 - alpha * q
    NB = 1
    while f ** NB * c >= c - h:
        NB += 1
    T_alpha = (81 * K) // (80 * alpha) + 1
    return {"NA": NA, "NB": NB, "U1": U1, "T": T_alpha}


def entry_and_itinerary(K, alpha, H, full_head=200):
    """Full exact vector until entry+full_head (closure asserted), then
    exact scalar map only (cheap) for the remainder. Returns (t_star,
    s_list of pre-update selections t=0..H-1)."""
    m, h, b, c, q, L = params(K)
    p = replay.exact_bin_probs(replay.near_cdf_exact, K)
    g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
    s_list = []
    t_star = None
    x = None
    t = 0
    while t < H:
        if t_star is None or t <= t_star + full_head:
            D = replay.cumsum_excess(p, K)
            s = replay.argmax_first(D)
            if t_star is not None:
                assert p[m - 1] == x, f"scalar closure failed K={K} t={t}"
            if t_star is None and s in (m - 1, m):
                t_star = t
                x = p[m - 1]
            if t_star is not None:
                x = ((1 - alpha) * x + alpha * L if x > h
                     else x + alpha * (c - x) * q)
            p = replay.update_state(p, g, s, alpha)
        else:
            s = m if x > h else m - 1
            x = ((1 - alpha) * x + alpha * L if x > h
                 else x + alpha * (c - x) * q)
        s_list.append(s)
        t += 1
    return t_star, s_list


def classify_two_test(s_list, K, alpha, H):
    m = params(K)[0]
    bd = bounds(K, alpha)
    tau0 = t_star_of(s_list, m) + max(bd["NB"], bd["U1"])
    ell_after = any(s == m - 1 for s in s_list[tau0 + 1:])
    return "CYCLE" if ell_after else "LOCK"


def t_star_of(s_list, m):
    return next(t for t, s in enumerate(s_list) if s in (m - 1, m))


def classify_frozen(s_list, K, H):
    m = params(K)[0]
    R_m = 0
    for s in reversed(s_list):
        if s == m:
            R_m += 1
        else:
            break
    t_star = t_star_of(s_list, m)
    post = set(s_list[t_star:])
    if R_m >= K:
        return "LOCK"
    if post == {m - 1, m}:
        return "CYCLE"
    return "OTHER"


def main():
    t0 = time.perf_counter()
    rows = []
    fails = []
    for alpha in CFG["alphas"]:
        # very small alpha windows are long (H ~ K/alpha); cap the lattice
        use_Ks = (CFG["Ks"] if alpha >= F(1, 50)
                  else CFG["Ks_cap_small_alpha"])
        for K in use_Ks:
            print(f"  [C] alpha={alpha} K={K}", flush=True)
            m, h, b, c, q, L = params(K)
            bd = bounds(K, alpha)
            H = bd["T"] + max(bd["NB"], bd["U1"]) + bd["NA"] + 2
            # cap for runtime: exact full vector to H is expensive for tiny
            # alpha; use scalar-only from entry with a short full-vector head
            Hp_pre = bd["T"] + max(bd["NB"], bd["U1"]) + \
                max(bd["NA"] + 2, K)
            t_star, s_list = entry_and_itinerary(K, alpha, max(H, Hp_pre))
            true_class = "LOCK" if K % 6 in (0, 1, 2) else "CYCLE"
            # C1: entry and dwell bounds actually checked on the window
            # (window = first H entries of the itinerary)
            win = s_list[:H]
            if t_star is None or t_star > bd["T"]:
                fails.append(("C1-entry", K, str(alpha), t_star, bd["T"]))
            import itertools
            runs = [(k, len(list(g))) for k, g in
                    itertools.groupby(win[t_star:])]
            # drop possibly truncated final run of the window
            if runs:
                runs = runs[:-1]
            a_runs = [ln for s_, ln in runs if s_ == m]
            b_runs = [ln for s_, ln in runs if s_ == m - 1]
            if a_runs:
                # first post-entry A-run (entry via m, x0 <= c) is bounded
                # by U1; all subsequent (post-landing) A-runs by NA alone
                first_is_A = bool(runs) and runs[0][0] == m
                later = a_runs[1:] if first_is_A else a_runs
                if first_is_A and a_runs[0] > bd["U1"]:
                    fails.append(("C1-A-first", K, str(alpha),
                                  a_runs[0], bd["U1"]))
                if later and max(later) > bd["NA"]:
                    fails.append(("C1-A-post", K, str(alpha),
                                  max(later), bd["NA"]))
            if b_runs and max(b_runs) > bd["NB"]:
                fails.append(("C1-B-dwell", K, str(alpha),
                              max(b_runs), bd["NB"]))
            full = True
            if full:
                # C4: LOCK never selects ell after t* + NB (within window)
                if true_class == "LOCK":
                    if any(s == m - 1 for s in win[t_star + bd["NB"] + 1:]):
                        fails.append(("C4", K, str(alpha)))
                # C2 two-test (within the declared window H)
                tau0 = t_star + max(bd["NB"], bd["U1"])
                ell_after = any(s == m - 1 for s in win[tau0 + 1:])
                two = "CYCLE" if ell_after else "LOCK"
                if two != true_class:
                    fails.append(("C2", K, str(alpha), two, true_class))
                # C3 frozen classifier: sufficient condition (NA < K =>
                # correct at H' = T + max(NB,U1) + max(NA+2, K)); when
                # NA >= K misfire is possible (observation only)
                Hp = bd["T"] + max(bd["NB"], bd["U1"]) + \
                    max(bd["NA"] + 2, K)
                fr = classify_frozen(s_list[:Hp], K, Hp)
                if fr != true_class and bd["NA"] < K:
                    fails.append(("C3", K, str(alpha), fr, true_class,
                                  bd["NA"]))
                rows.append({"alpha": str(alpha), "K": K, "H": H, "Hp": Hp,
                             "t_star": t_star, "NA": bd["NA"],
                             "NB": bd["NB"], "U1": bd["U1"], "T": bd["T"],
                             "two_test": two, "frozen": fr,
                             "true": true_class})
    # counterexample reproduction (frozen, exact)
    ce = CFG["counterexample"]
    m = params(ce["K"])[0]
    t_star, s_list = entry_and_itinerary(ce["K"], ce["alpha"], ce["H"])
    fr = classify_frozen(s_list, ce["K"], ce["H"])
    bd = bounds(ce["K"], ce["alpha"])
    tau0 = t_star + max(bd["NB"], bd["U1"])
    two = "CYCLE" if any(s == m - 1 for s in s_list[tau0 + 1:]) else "LOCK"
    ce_res = {"frozen": fr, "two_test": two, "NA": bd["NA"],
              "expected_frozen": "LOCK", "expected_two": "CYCLE"}
    ce_ok = (fr == "LOCK" and two == "CYCLE" and bd["NA"] >= ce["K"])

    res = {"rows": len(rows), "fails": fails, "counterexample": ce_res,
           "counterexample_ok": ce_ok,
           "elapsed_s": round(time.perf_counter() - t0, 1)}
    with (OUT / "c_rows.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    res["git"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                                capture_output=True, text=True).stdout.strip()
    res["output_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(OUT.glob("*.csv"))}
    (OUT / "c_summary.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8")
    print(json.dumps(res, indent=2, ensure_ascii=False, default=str))
    return 0 if (not fails and ce_ok) else 1


if __name__ == "__main__":
    raise SystemExit(main())
