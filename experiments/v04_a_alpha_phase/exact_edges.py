# -*- coding: utf-8 -*-
"""v04-A exact edge certification: window-exact rechecks of the two common-
window edges found by the float scan, plus dyadic-rational bisection of the
critical alphas.  Window evidence (T=100 exact replay), NOT yet all-time
certificates; residence flag recorded for later tail upgrade.

Edges probed:
  LEFT  : K=10, alpha_c(10) in (59/240, 119/480), binding violation t=3
  RIGHT : K=6 , alpha_c^+(6) in (17/60, 7/24),     binding violation t=1
  RIGHT2: K=12, K=13 edges in (3/10, 13/40)        (binding t=1)
"""
from __future__ import annotations

import json
from fractions import Fraction as F
from pathlib import Path

T = 100


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


def replay(K: int, alpha: F, T: int = T):
    p0 = [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]
    g = [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]
    m = -(-5 * K // 6)
    states = [p0]
    sels = []
    for _ in range(T):
        p = states[-1]
        acc, best = F(0), 0
        Dbest = None
        for j, pj in enumerate(p):
            acc += pj
            d = acc - F(j + 1, K)
            if Dbest is None or d > Dbest:
                Dbest, best = d, j
        s = best + 1
        M = alpha * sum(p[:s])
        states.append([((1 - alpha) * pj if j < s else pj) + M * gj
                       for j, (pj, gj) in enumerate(zip(p, g))])
        sels.append(s)
    U = [K * sum((y - F(1, K)) ** 2 for y in st) for st in states]
    return U, sels, m


def gap_profile(K: int, alpha: F):
    U, sels, m = replay(K, alpha)
    t_star = next((t for t, s in enumerate(sels) if s in (m - 1, m)), None)
    residence = t_star is not None and all(s in (m - 1, m) for s in sels[t_star:])
    min_t3 = min(range(3, T + 1), key=lambda t: U[t])
    return {
        "U1_minus_U2": U[1] - U[2],
        "U3_minus_U2": U[3] - U[2],
        "min_tail_minus_U2": U[min_t3] - U[2],
        "argmin_tail": min_t3,
        "t_star": t_star, "residence": residence,
        "sels_prefix": sels[:5],
    }


def fmt(x: F) -> str:
    f = float(x)
    return f"{x.numerator}/{x.denominator} ({f:.3e})"


def bisect(K: int, lo: F, hi: F, side: str, iters: int = 18):
    """side='left': violation (tail min < U2) for alpha <= edge
       side='right': violation (U1 < U2) for alpha >= edge"""
    def bad(a: F) -> bool:
        p = gap_profile(K, a)
        return (p["min_tail_minus_U2"] < 0 if side == "left"
                else p["U1_minus_U2"] < 0)
    if side == "left":
        assert bad(lo) and not bad(hi), "bracket invalid"
    else:
        assert (not bad(lo)) and bad(hi), "bracket invalid"
    for _ in range(iters):
        mid = (lo + hi) / 2
        if side == "left":
            # bad for alpha <= edge: bad(mid) -> edge is above mid
            if bad(mid):
                lo = mid
            else:
                hi = mid
        else:
            # bad for alpha >= edge: bad(mid) -> edge is below mid
            if bad(mid):
                hi = mid
            else:
                lo = mid
    return lo, hi


def main():
    out = {}
    print("== exact rechecks (T=100 exact replay) ==")
    for K, a in [(10, F(59, 240)), (10, F(119, 480)),
                 (6, F(7, 24)), (6, F(17, 60))]:
        p = gap_profile(K, a)
        tag = ("VIOLATION" if (p["min_tail_minus_U2"] < 0
                               or p["U1_minus_U2"] < 0) else "clean")
        print(f"K={K} a={a} -> {tag}: U1-U2={fmt(p['U1_minus_U2'])} "
              f"U3-U2={fmt(p['U3_minus_U2'])} tailmin-U2="
              f"{fmt(p['min_tail_minus_U2'])}@t={p['argmin_tail']} "
              f"t*={p['t_star']} res={p['residence']} sels={p['sels_prefix']}")
        out[f"recheck_{K}_{a.numerator}_{a.denominator}"] = {
            k: (str(v) if not isinstance(v, list) else v)
            for k, v in p.items()}

    print("\n== dyadic bisection of critical alphas ==")
    for K, lo, hi, side, name in [
            (10, F(59, 240), F(119, 480), "left", "alpha_c(10) left edge"),
            (6, F(17, 60), F(7, 24), "right", "alpha_c+(6) right edge"),
            (12, F(3, 10), F(13, 40), "right", "alpha_c+(12)"),
            (13, F(3, 10), F(13, 40), "right", "alpha_c+(13)"),
            (9, F(23, 96), F(29, 120), "left", "alpha_c(9) left edge")]:
        a, b = bisect(K, lo, hi, side)
        val = (a + b) / 2
        print(f"{name}: in ({a}, {b})  ~ {float(val):.9f}")
        out[name] = {"lo": str(a), "hi": str(b), "mid_float": float(val)}

    dest = Path(__file__).with_name("exact_edges.json")
    dest.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", dest)


if __name__ == "__main__":
    main()
