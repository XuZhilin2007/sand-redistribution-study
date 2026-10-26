# -*- coding: utf-8 -*-
"""v04-C1: minimal tail-leakage probe, EXACT arithmetic (R2-repaired).

Model change vs canonical v0.3.0 baseline (single change):
    respray kernel g_eps = (1-eps) * g_XO + eps * u,  u_i = 1/K.
Init canonical, alpha = 1/4, smallest-index argmax, exact Fractions.

NATURAL STOPPING (review R2): before each update, compute
d = max_j D_j; if d <= tau = 1e-12 the process STOPS: state frozen,
no further active updates.  Selections at the stop moment (e.g. argmax
at j=K with d=0) are stop events, not legal active beyond-m updates.

R5 machinery (machine-verified identity): while all active selections
so far satisfy s <= m, for j >= m:
    D_{t,j} = D_{0,j} - eps*(1-y_j)*R_t,   R_t = sum_{k<t} M_k.
Hence a first legal active selection s>m requires (y-b)[eps*R_t -
(b+y-1)] > 0 with y>b, i.e. eps*R_t > 2b+h-1; since R_t <= alpha*t,
NECESSARY bound: t > (2b+h-1)/(alpha*eps).  The script asserts the
identity exactly on every step and checks the bound at the first
beyond-m event.

Diagnostics per (K, eps): stop time & d at stop, active steps, U
series to stop, argmin over legal states, loss of stopping at 2,
first legal beyond-m selection, exact necessary-bound check.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path

TAU = F(1, 10 ** 12)
T = 120
CASES = [(K, e) for K in (6, 9, 12) for e in
         (F(0), F(1, 100), F(1, 20), F(1, 10), F(1, 4))]
BOUND_CHECKS = {(12, F(1, 10))}   # review's probe point


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


def run(K: int, eps: F, T: int = T):
    p0 = [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]
    g_xo = [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]
    g = [(1 - eps) * gj + eps * F(1, K) for gj in g_xo]
    m = -(-5 * K // 6)
    b = F(m, K)
    h = F(1, K)
    c = 2 * b - b * b
    alpha = F(1, 4)
    states = [p0]
    sels = []
    R = F(0)                    # cumulative removed mass (legal steps)
    identity_ok = True
    bound_ok = None             # filled at first beyond-m event
    stopped_at = None
    d_stop = None
    for t in range(T):
        p = states[-1]
        acc, Dbest, best = F(0), None, 0
        D = []
        for j, pj in enumerate(p):
            acc += pj
            d = acc - F(j + 1, K)
            D.append(d)
            if Dbest is None or d > Dbest:
                Dbest, best = d, j
        # natural stop BEFORE updating
        if Dbest <= TAU:
            stopped_at, d_stop = t, Dbest
            break
        s = best + 1
        M = alpha * sum(p[:s])
        # R5 exact identity check (holds while all past s <= m)
        if all(x <= m for x in sels):
            y_list = [F(j + 1, K) for j in range(m - 1, K)]
            base = [yy * (1 - yy) for yy in y_list]   # canonical D_0 = y(1-y)
            for idx, j in enumerate(range(m - 1, K)):
                lhs = sum(p[:j + 1]) - F(j + 1, K)
                if lhs != base[idx] - eps * (1 - y_list[idx]) * R:
                    identity_ok = False
        if s > m and bound_ok is None:
            # first legal beyond-m: necessary bound t > (2b+h-1)/(a*eps)
            need = (2 * b + h - 1) / (alpha * eps)
            bound_ok = (t > need, t, f"{need.numerator}/{need.denominator}")
        sels.append(s)
        states.append([((1 - alpha) * pj if j < s else pj) + M * gj
                       for j, (pj, gj) in enumerate(zip(p, g))])
        R += M
    U = [K * sum((y - F(1, K)) ** 2 for y in st) for st in states]
    t_arg = min(range(len(states)), key=lambda t: U[t])
    loss2 = (U[2] - U[t_arg]) / U[2] if len(states) > 2 else F(0)
    drift = [sum(st[:m]) - c for st in states]
    return {
        "K": K, "eps": f"{eps.numerator}/{eps.denominator}",
        "m": m, "steps_recorded": len(states) - 1,
        "stopped_at": stopped_at,
        "d_stop": None if d_stop is None else
        f"{d_stop.numerator}/{d_stop.denominator}",
        "identity_ok": identity_ok,
        "first_beyond_m": next((s for s in sels if s > m), None),
        "bound_check": bound_ok,
        "argmin": t_arg, "loss2_float": float(loss2),
        "drift_Fm_last": float(drift[-1]),
        "U_last": f"{U[-1].numerator}/{U[-1].denominator}",
        "sels_prefix": sels[:12],
    }


def main() -> int:
    out = []
    for K, e in CASES:
        r = run(K, e)
        out.append(r)
        bc = r["bound_check"]
        print(f"K={K:>2} eps={r['eps']:>6} stop@{r['stopped_at']} "
              f"steps={r['steps_recorded']} argmin={r['argmin']} "
              f"loss2={r['loss2_float']:+.6f} beyond_m={r['first_beyond_m']} "
              f"id_ok={r['identity_ok']} drift={r['drift_Fm_last']:+.2e}")
        if bc is not None:
            passed, t, need = bc
            print(f"    R5 bound at first beyond-m: t={t} > {need}"
                  f" required -> {passed}")
    hard = [r for r in out if not r["identity_ok"]] + \
           [r for r in out if r["bound_check"] and not r["bound_check"][0]]
    dest = Path(__file__).with_name("leak_probe.json")
    dest.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", dest)
    print("HARD FAILURES:", len(hard))
    return 1 if hard else 0


if __name__ == "__main__":
    sys.exit(main())
