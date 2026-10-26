# -*- coding: utf-8 -*-
"""v04-A: exact piecewise analysis of U1-U2 for K=6 in alpha.

D_{1,j}(alpha) is affine in alpha (single update of the fixed canonical
D_0 with selection s_0 fixed), so every s_1 switch point is a RATIONAL
alpha.  On each piece s_1 is constant and U1-U2 is a polynomial in alpha
with rational coefficients (computed exactly with a mini poly engine).
This locates the right edge alpha_c^+(6) exactly: either a rational
switch point or an algebraic root inside a piece.
"""
from __future__ import annotations

from fractions import Fraction as F

K = 6


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


# --- tiny polynomial engine over Fractions -----------------------------
def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(n)]


def pmul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def peval(a, x):
    r = F(0)
    for c in reversed(a):
        r = r * x + c
    return r


# --- model data ---------------------------------------------------------
g = [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]
h = F(1, K)
p0 = [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]
D0 = []
acc = F(0)
for j, pj in enumerate(p0):
    acc += pj
    D0.append(acc - F(j + 1, K))
# s0 = smallest argmax of D0
s0 = max(range(K), key=lambda j: (D0[j], -j)) + 1
Fs0 = sum(p0[:s0])
print(f"K=6  s0={s0}  Fs0={Fs0}")

# D_{1,j}(alpha) affine: slope/intercept
sl = []
ic = []
for j in range(K):
    if j + 1 <= s0:
        # D' = D + alpha*(Fs0*Gj - jh - D)
        Gj = sum(g[:j + 1])
        base = D0[j]
        sl.append(Fs0 * Gj - F(j + 1, K) - base)
        ic.append(base)
    else:
        Gj = sum(g[:j + 1])
        sl.append(-Fs0 * (1 - Gj))
        ic.append(D0[j])
# s1 selection as function of alpha: argmax_j ic_j + alpha*sl_j (smallest idx)
switches = []
for a in range(K):
    for b in range(a + 1, K):
        if sl[a] != sl[b]:
            x = (ic[b] - ic[a]) / (sl[a] - sl[b])
            if F(1, 5) < x < F(9, 10):
                switches.append(x)
switches = sorted(set(switches))
print("s1 switch points (rational) in (0.2, 0.9):")
for x in switches:
    # which pair crosses and who wins where
    print(f"  {x} = {float(x):.9f}")


def s1_at(a: F) -> int:
    vals = [ic[j] + a * sl[j] for j in range(K)]
    best = 0
    for j in range(1, K):
        if vals[j] > vals[best]:
            best = j
    return best + 1


# --- U1, U2 as polynomials on a fixed (s0, s1) piece --------------------
def step_poly(p, s, alpha_poly):
    """p: list of polys in alpha (bin masses). Returns stepped masses and
    M as poly. alpha_poly = [0, 1]."""
    Fs = [F(0)]
    for i in range(s):
        Fs = padd(Fs, p[i])
    M = pmul(alpha_poly, Fs)
    out = []
    for j in range(K):
        base = pmul([F(1)], p[j]) if j + 1 > s else padd(
            p[j], pmul([F(0), F(-1)], p[j]))
        out.append(padd(base, pmul(M, [g[j]])))
    return out


def U_of(p):
    """U as poly: K * sum (p_i - h)^2."""
    s = [F(0)]
    for i in range(K):
        d = padd(p[i], [-h])
        s = padd(s, pmul(d, d))
    return [K * c for c in s]


alpha = [F(0), F(1)]
p0p = [[v] for v in p0]
p1p = step_poly(p0p, s0, alpha)
p2p = step_poly(p1p, s1_at(F(1, 4)), alpha)   # placeholder; recomputed below

print("\npieces (alpha interval, s1), U1-U2 sign at ends:")
pts = [F(1, 5)] + switches + [F(9, 10)]
for lo, hi in zip(pts[:-1], pts[1:]):
    if hi <= F(1, 4) or lo >= F(2, 5):
        continue  # focus on the edge region (0.25, 0.4)
    mid = (lo + hi) / 2
    s1 = s1_at(mid)
    p1 = step_poly(p0p, s0, alpha)
    p2 = step_poly(p1, s1, alpha)
    gU1 = U_of(p1)
    gU2 = U_of(p2)
    diff = padd(gU1, [-c for c in gU2])
    v_lo, v_hi, v_mid = peval(diff, lo), peval(diff, hi), peval(diff, mid)
    print(f"  a in ({float(lo):.6f}, {float(hi):.6f}) s1={s1}: "
          f"U1-U2 = {float(v_lo):+.3e} .. {float(v_mid):+.3e} .. "
          f"{float(v_hi):+.3e}   deg={len(diff)-1}")
    if (v_lo > 0 > v_hi) or (v_lo < 0 < v_hi):
        a, b = lo, hi
        for _ in range(60):
            m = (a + b) / 2
            if (peval(diff, a) > 0) == (peval(diff, m) > 0):
                a = m
            else:
                b = m
        print(f"    ROOT at alpha ~ {(a+b)/2} = {float((a+b)/2):.12f}")
