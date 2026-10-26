# -*- coding: utf-8 -*-
"""v04-A anchor replay: 25 exact (K, alpha) anchors, 0..60 rounds, exact
Fraction arithmetic, self-contained model (mirrors the S3 definitions; the
step() semantics follow the v0.3.0 independent verifier conventions).

Model (only change vs canonical v0.3.0 baseline: general alpha):
  p0_i = (2K-2i+1)/K^2                      (canonical init)
  g_i  = G(i/K)-G((i-1)/K), G = frozen XO CDF
  s_t  = smallest index maximizing D_{t,j} = F_{t,j} - j/K
  p_{t+1,i} = (1-alpha) p_{t,i} + alpha*F_{t,s}*g_i   (i <= s)
              p_{t,i}         + alpha*F_{t,s}*g_i     (i >  s)
  U(p) = K * sum_i (p_i - 1/K)^2

Outputs per (K, alpha): exact U_t series 0..60, earliest window argmin,
t_* = first boundary selection (s_t in {m-1, m}), whether U(t) < U(2)
anywhere in the window, and the K=9 alpha=1/5 counterexample values.
Exit 0 always (this is a probe, not a gate); JSON carries the evidence.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path

T_WINDOW = 60
KS = [9, 10, 11, 12, 65]
ALPHAS = [F(1, 8), F(1, 5), F(1, 4), F(1, 3), F(1, 2)]


def G(x: F) -> F:
    if x <= 0:
        return F(0)
    if x < F(1, 3):
        return x / 2
    if x <= F(2, 3):
        return 3 * x / 2 - F(1, 3)
    if x <= F(5, 6):
        return 2 * x - F(2, 3)
    return F(1)


def bins_p0(K: int):
    return [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]


def bins_g(K: int):
    return [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]


def step(p, g, K, alpha):
    D = []
    acc = F(0)
    best = 0
    for j, pj in enumerate(p):
        acc += pj
        d = acc - F(j + 1, K)
        D.append(d)
        if d > D[best]:
            best = j
    s = best + 1
    F_s = F(1, 1) * sum(p[:s])
    M = alpha * F_s
    out = []
    for j, (pj, gj) in enumerate(zip(p, g)):
        base = (1 - alpha) * pj if j < s else pj
        out.append(base + M * gj)
    return out, s


def run_anchor(K: int, alpha: F):
    p0 = bins_p0(K)
    g = bins_g(K)
    m = -(-5 * K // 6)
    states = [p0]
    sels = []
    for _ in range(T_WINDOW):
        p2, s = step(states[-1], g, K, alpha)
        sels.append(s)
        states.append(p2)
    U = [K * sum((y - F(1, K)) ** 2 for y in st) for st in states]
    # earliest window argmin
    t_arg = min(range(T_WINDOW + 1), key=lambda t: (U[t], t))
    # first boundary selection time (selection indices are 1-based)
    t_star = next((t for t, s in enumerate(sels) if s in (m - 1, m)), None)
    viol = [t for t in range(T_WINDOW + 1) if U[t] < U[2]]
    # selection word prefix (until it stabilises into {m-1, m} forever)
    return {
        "K": K,
        "alpha": f"{alpha.numerator}/{alpha.denominator}",
        "m": m,
        "t_star": t_star,
        "t_argmin_window": t_arg,
        "U_argmin_over_U2": f"{(U[t_arg] / U[2]).numerator}/"
                            f"{(U[t_arg] / U[2]).denominator}",
        "violations_of_U2": viol,
        "sel_word_prefix": sels[:12],
        "U_series": [f"{u.numerator}/{u.denominator}" for u in U],
    }


def main() -> int:
    out = []
    for K in KS:
        for alpha in ALPHAS:
            rec = run_anchor(K, alpha)
            out.append(rec)
            print(f"K={K:>3} a={rec['alpha']:>4} t*={rec['t_star']} "
                  f"argmin={rec['t_argmin_window']} "
                  f"viol(U<U2)={rec['violations_of_U2'][:4]} "
                  f"word={rec['sel_word_prefix'][:8]}")
    # highlight the precise counterexample K=9, alpha=1/5
    ce = next(r for r in out if r["K"] == 9 and r["alpha"] == "1/5")
    U2 = F(*map(int, ce["U_series"][2].split("/")))
    U3 = F(*map(int, ce["U_series"][3].split("/")))
    print("\n[K=9, alpha=1/5] U(p_2) =", U2, " U(p_3) =", U3)
    print("  U(p_3) - U(p_2) =", U3 - U2, " -> negative:", U3 < U2)
    dest = Path(__file__).with_name("anchors25_exact.json")
    dest.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
