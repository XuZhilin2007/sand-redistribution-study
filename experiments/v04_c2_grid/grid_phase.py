# -*- coding: utf-8 -*-
"""v04-C2: grid-phase probe (exact arithmetic).

Continuous model fixed: init CDF A(x)=2x-x^2, respray CDF G = XO.
Grid 0=x0<...<xK=1 given by warp phi(x)=(1-lam)x+lam*x^2 (increasing,
phi(0)=0, phi(1)=1).  Cell masses p0_i = A(x_i)-A(x_{i-1}); target
u_i = x_i - x_{i-1} (cell LENGTH, not 1/K); kernel g_i = G(x_i)-G(x_{i-1}).
Metric U_w = sum (p_i-u_i)^2 / u_i  (weighted; equals canonical U on
uniform grids where u_i = 1/K).

Discrepancy D_j = F_j - x_j (cumulative mass minus cumulative target =
cumulative LENGTH).  Selection = smallest argmax.  Removal alpha=1/4 of
selected prefix, respray M = alpha*F_s over g.

Checks and observables:
  (0) consistency: nonneg p, sums = 1, uniform-grid reduction (lam=0)
      reproduces canonical m/q/conservation F_{t,m}=c.
  (1) uniform lam=0, K=12..17: the mod-6 phase r -> boundary-cell phase
      sigma = (5/6 - x_{m-1})/h_cell, word type, two-round optimality.
  (2) warped lam in {1/10, 1/4}, K=12..17: does the 6-phase structure
      persist (same word classes) or dissolve? sigma moves continuously.
  (3) K=12 fixed, lam swept finely: word switches as sigma crosses cell?
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path

T = 100
ALPHA = F(1, 4)


def Gc(x: F) -> F:
    if x <= 0:
        return F(0)
    if x < F(1, 3):
        return x / 2
    if x <= F(2, 3):
        return 3 * x / 2 - F(1, 3)
    if x < F(5, 6):
        return 2 * x - F(2, 3)
    return F(1)


def Ac(x: F) -> F:
    return 2 * x - x * x


def grid(K: int, lam: F):
    xs = []
    for i in range(K + 1):
        x = F(i, K)
        xs.append((1 - lam) * x + lam * x * x)
    return xs


def run(K: int, lam: F, T: int = T):
    xs = grid(K, lam)
    p = [Ac(xs[i]) - Ac(xs[i - 1]) for i in range(1, K + 1)]
    u = [xs[i] - xs[i - 1] for i in range(1, K + 1)]
    g = [Gc(xs[i]) - Gc(xs[i - 1]) for i in range(1, K + 1)]
    assert all(v >= 0 for v in p) and all(v > 0 for v in u)
    assert sum(p) == 1 and sum(u) == 1 and sum(g) == 1
    # support end cell: smallest m with G(x_m)=1
    m = next(j for j in range(1, K + 1) if Gc(xs[j]) == 1)
    sigma = (F(5, 6) - xs[m - 1]) / (xs[m] - xs[m - 1])
    states = [p]
    sels = []
    stopped_at = None
    for t in range(T):
        q = states[-1]
        acc, Dbest, best = F(0), None, 0
        for j, pj in enumerate(q):
            acc += pj
            d = acc - xs[j + 1]
            if Dbest is None or d > Dbest:
                Dbest, best = d, j
        if Dbest <= F(1, 10 ** 12):        # natural stop BEFORE updating
            stopped_at = t
            break
        s = best + 1
        M = ALPHA * sum(q[:s])
        states.append([((1 - ALPHA) * pj if j < s else pj) + M * gj
                       for j, (pj, gj) in enumerate(zip(q, g))])
        sels.append(s)
    Uw = [sum((st[i] - u[i]) ** 2 / u[i] for i in range(K))
          for st in states]
    t_arg = min(range(T + 1), key=lambda t: Uw[t])
    word = sels[-24:]
    # reduce word to its minimal period if periodic
    per = next((pp for pp in range(1, 13)
                if all(word[i] == word[i + pp]
                       for i in range(len(word) - pp))), None)
    t_star = next((t for t, s in enumerate(sels) if s in (m - 1, m)), None)
    residence = t_star is not None and all(
        s in (m - 1, m) for s in sels[t_star:])
    beyond = [s for s in sels if s > m]
    # conservation check of F at m
    Fm = [sum(st[:m]) for st in states]
    return {
        "K": K, "lam": f"{lam.numerator}/{lam.denominator}",
        "m": m, "sigma": f"{sigma.numerator}/{sigma.denominator}",
        "sigma_f": float(sigma), "t_star": t_star, "residence": residence,
        "stopped_at": stopped_at,
        "beyond_m_count": len(beyond),
        "Fm_conserved": Fm[T] == Fm[0],
        "argmin_Uw": t_arg,
        "Uw2_is_min": t_arg == 2,
        "word_period": per, "word_tail": word[:8],
        "Uw_series_head": [f"{v.numerator}/{v.denominator}" for v in Uw[:6]],
    }


def main() -> int:
    out = []
    # (0) uniform reduction sanity: K=12 lam=0 must give canonical facts
    r = run(12, F(0))
    print(f"[sanity] K=12 lam=0: m={r['m']} (canonical m=10), "
          f"Fm conserved={r['Fm_conserved']}, argmin={r['argmin_Uw']} "
          f"(canonical t_off=2), word={r['word_tail']} per={r['word_period']}")
    assert r["m"] == 10 and r["Fm_conserved"] and r["argmin_Uw"] == 2
    print("  -> uniform-grid reduction consistent with canonical\n")

    print("(1) uniform lam=0, K=12..17 (mod-6 phases):")
    for K in range(12, 18):
        r = run(K, F(0))
        out.append(r)
        print(f"  K={K} m={r['m']:>2} sigma={r['sigma_f']:.4f} "
              f"res={r['residence']} argmin={r['argmin_Uw']} "
              f"per={r['word_period']} word={r['word_tail']}")
    print("\n(2) warped lam=1/10 and 1/4, K=12..17:")
    for lam in (F(1, 10), F(1, 4)):
        for K in range(12, 18):
            r = run(K, lam)
            out.append(r)
            print(f"  lam={float(lam):.2f} K={K} m={r['m']:>2} "
                  f"sigma={r['sigma_f']:.4f} res={r['residence']} "
                  f"argmin={r['argmin_Uw']} per={r['word_period']} "
                  f"word={r['word_tail']} beyond={r['beyond_m_count']}")
    print("\n(3) K=12, lam swept:")
    for lam in [F(i, 80) for i in range(0, 41, 2)]:
        r = run(12, lam)
        out.append(r)
        print(f"  lam={float(lam):.4f} m={r['m']} sigma={r['sigma_f']:.4f} "
              f"res={r['residence']} argmin={r['argmin_Uw']} "
              f"per={r['word_period']} word={r['word_tail'][:6]}")
    dest = Path(__file__).with_name("grid_phase.json")
    dest.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
