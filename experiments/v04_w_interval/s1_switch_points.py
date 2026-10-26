# -*- coding: utf-8 -*-
"""v04 phase-2 companion: exact s1 switch points near 1/4 for the
W1-OPEN Ks (those whose D1-line switch lies within the radius ladder).

D_1,j(alpha) is affine; every pairwise crossing is a RATIONAL alpha.
We list all crossings in (1/4 - 1/125, 1/4 + 1/125), identify which
involve the running argmax (real switches), and report the exact
distance |alpha* - 1/4| -- the per-K certification radius obstacle.
"""
from __future__ import annotations

import sys
from fractions import Fraction as F

sys.path.insert(0, ".")


def Gc(x: F) -> F:
    if x <= 0:
        return F(0)
    if x < F(1, 3):
        return x / 2
    if x <= F(2, 3):
        return 3 * x / 2 - F(1, 3)
    if x <= F(5, 6):
        return 2 * x - F(2, 3)
    return F(1)


def lines(K: int):
    p0 = [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]
    g = [Gc(F(j, K)) - Gc(F(j - 1, K)) for j in range(1, K + 1)]
    D0, acc = [], F(0)
    for j, pj in enumerate(p0):
        acc += pj
        D0.append(acc - F(j + 1, K))
    s0 = 1 + max(range(K), key=lambda j: (D0[j], -j))
    Fs0 = sum(p0[:s0])
    ic, sl = [], []
    for j in range(K):
        Gj = sum(g[:j + 1])
        ic.append(D0[j])
        sl.append((Fs0 * Gj - F(j + 1, K) - D0[j]) if j + 1 <= s0
                  else -Fs0 * (1 - Gj))
    return ic, sl


def real_switches(K: int, halfwidth=F(1, 125)):
    ic, sl = lines(K)
    a0 = F(1, 4)
    # running argmax at a0 and at interval ends
    def arg(a):
        return max(range(K), key=lambda j: (ic[j] + a * sl[j], -j))
    switches = []
    for a_end in (a0 - halfwidth, a0 + halfwidth):
        cur = arg(a0)
        # walk is not needed: collect ALL crossings, keep those where
        # the two lines are the top-two locally (verified by sampling)
        for j in range(K):
            for jj in range(j + 1, K):
                if sl[j] == sl[jj]:
                    continue
                x = (ic[jj] - ic[j]) / (sl[j] - sl[jj])
                if a0 - halfwidth < x < a0 + halfwidth:
                    # real switch iff at x-eps the argmax is one of them
                    # and at x+eps it is the other: check argmax identity
                    lo_arg = arg(x - F(1, 10 ** 9))
                    hi_arg = arg(x + F(1, 10 ** 9))
                    if {lo_arg, hi_arg} == {j, jj}:
                        switches.append((x, j + 1, jj + 1))
    switches = sorted(set(switches))
    return switches, arg(a0) + 1


def main():
    for K in [16, 32, 96, 107, 110, 124, 129, 132, 143, 146, 157, 160,
              171, 174, 179, 182, 185, 188, 193, 196, 199, 202]:
        sw, s1 = real_switches(K)
        if not sw:
            print(f"K={K}: no REAL switch inside +-1/125 (s1={s1})")
            continue
        dists = [(abs(x - F(1, 4)), x, a, b) for x, a, b in sw]
        dmin, x, a, b = min(dists)
        print(f"K={K:>3}: {len(sw)} real s1 switch(es) in window; "
              f"nearest |a*-1/4| = {dmin} = {float(dmin):.3e} "
              f"(a*={float(x):.9f}, lines {a}<->{b}); "
              f"s1(1/4)={s1}")


if __name__ == "__main__":
    main()
