# -*- coding: utf-8 -*-
"""v04-D: non-canonical initial states, exact probe (alpha=1/4, XO).

R2/R3-repaired: natural stopping (d<=tau freezes the process BEFORE any
update) is honoured, and the P3 family is a genuinely nonzero
boundary-preserving perturbation (prefix-internal move bin m-3 -> m-1,
which keeps F_m = c exactly); the actual L1 change is recorded and a
zero change at delta>0 is a hard failure.

Families (all legal: nonneg, sum 1):
  P1 bump        : move delta from bin K to bin 1 (changes F_m).
  P2 tilt        : p += delta*(2i-K+1)/K (sum preserved).
  P3 bfix (new)  : move delta from bin m-3 to bin m-1 (prefix-internal;
                   F_m, c, a and the whole suffix untouched).
Observables: d_0 (natural-stop check at t=0), stopped_at, t*, window
argmin of U over LEGAL states, U(t)>=U(2) violations, L1 perturbation.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path

T = 60
TAU = F(1, 10 ** 12)


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


def traj(p0, K, T=T, alpha=F(1, 4)):
    g = [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]
    m = -(-5 * K // 6)
    states = [p0]
    sels = []
    stopped_at = None
    for t in range(T):
        p = states[-1]
        acc, Dbest, best = F(0), None, 0
        for j, pj in enumerate(p):
            acc += pj
            d = acc - F(j + 1, K)
            if Dbest is None or d > Dbest:
                Dbest, best = d, j
        if Dbest <= TAU:            # natural stop BEFORE updating
            stopped_at = t
            break
        s = best + 1
        M = alpha * sum(p[:s])
        states.append([((1 - alpha) * pj if j < s else pj) + M * gj
                       for j, (pj, gj) in enumerate(zip(p, g))])
        sels.append(s)
    return states, sels, m, stopped_at


def d0_of(p0, K):
    acc, dmax = F(0), None
    for j, pj in enumerate(p0):
        acc += pj
        d = acc - F(j + 1, K)
        if dmax is None or d > dmax:
            dmax = d
    return dmax


def analyze(p0, K, ref_p0=None):
    l1 = (sum(abs(a - b) for a, b in zip(p0, ref_p0))
          if ref_p0 is not None else F(0))
    states, sels, m, stopped_at = traj(p0, K)
    U = [K * sum((y - F(1, K)) ** 2 for y in st) for st in states]
    if len(states) < 3:
        # natural stop before two updates: no legal U(2) exists; the
        # process freezes at (or before) the initial state
        return {"stopped_at": stopped_at,
                "d0": str(d0_of(p0, K)), "d0_pos": d0_of(p0, K) > TAU,
                "L1_vs_canonical": str(l1),
                "argmin": 0, "viol": [], "n_viol": 0,
                "t_star": None, "residence": None,
                "note": "natural stop before t=2 (d<=tau); "
                        "no-intervention is forced, not chosen"}
    t_arg = min(range(len(states)), key=lambda t: U[t])
    viol = [t for t in range(len(states)) if U[t] < U[2]]
    t_star = next((t for t, s in enumerate(sels) if s in (m - 1, m)), None)
    res = t_star is not None and all(
        s in (m - 1, m) for s in sels[t_star:])
    return {"stopped_at": stopped_at,
            "d0": str(d0_of(p0, K)), "d0_pos": d0_of(p0, K) > TAU,
            "L1_vs_canonical": str(l1),
            "argmin": t_arg, "viol": viol[:3], "n_viol": len(viol),
            "t_star": t_star, "residence": res}


def legal(p):
    return all(v >= 0 for v in p) and sum(p) == 1


def main() -> int:
    out = []
    hard = []
    for K in (9, 12):
        p0 = [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]
        m = -(-5 * K // 6)
        print(f"== K={K} (m={m}) ==")
        r = analyze(p0, K)
        print(f"  canonical: {r}")
        out.append({"K": K, "family": "canonical", "delta": "0", **r})
        for delta in (F(1, 100), F(1, 20), F(1, 10)):
            # P1 bump bin1 <- binK
            p = list(p0)
            p[0] += delta
            p[K - 1] -= delta
            if legal(p):
                r = analyze(p, K, p0)
                out.append({"K": K, "family": "P1_bump", "delta": str(delta),
                            **r})
                print(f"  P1 d={delta}: {r}")
            # P2 tilt
            p = [p0[i] + delta * F(2 * i - K + 1, K) for i in range(K)]
            if legal(p):
                r = analyze(p, K, p0)
                out.append({"K": K, "family": "P2_tilt", "delta": str(delta),
                            **r})
                print(f"  P2 d={delta}: {r}")
            # P3 boundary-fixed: prefix-internal move bin m-3 -> bin m-1
            # (0-based indices m-4 -> m-2); requires m >= 4 and enough mass
            if m - 4 >= 0 and p0[m - 4] >= delta:
                p = list(p0)
                p[m - 4] -= delta
                p[m - 2] += delta
                assert legal(p) and sum(p[:m]) == sum(p0[:m])
                r = analyze(p, K, p0)
                if delta > 0 and r["L1_vs_canonical"] == "0":
                    hard.append(("P3-zero-change", K, str(delta)))
                out.append({"K": K, "family": "P3_bfix", "delta": str(delta),
                            **r})
                print(f"  P3 d={delta}: {r}")
            else:
                print(f"  P3 d={delta}: SKIPPED (mass/index constraint)")
    dest = Path(__file__).with_name("init_probe.json")
    dest.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", dest)
    print("HARD FAILURES:", hard)
    return 1 if hard else 0


if __name__ == "__main__":
    sys.exit(main())
