# -*- coding: utf-8 -*-
"""INDEPENDENT second implementation for Theorem E (CYCLE all-K all-time
optimal stopping), owner review protocol; review-10 repaired: the
quadratic identity carries a factor 2 on the cross term,
    U(p_t) = Q(x_t) + 2 z Phi(x_t) + K E2 z^2,  z = (3/4)^{t-2},
and this verifier now DIRECTLY back-substitutes the assembled quadratic
against the true U at t = 2, 3, 4 on every sampled K (the missing check
that let the factor-2 slip through in the previous round).  It also
checks the square-completion sufficient conditions S1-S4 and the
THEOREM STATEMENT itself on trajectories.

Imports NOTHING from run_e.py, sand_m0, or the frozen replay module --
rebuilds everything from the S3 model definitions with its OWN exact
engine.  Exit code 0 = all sampled K verified; nonzero = a real
counterexample or failure (printed with exact values).
"""
from __future__ import annotations

import sys
from fractions import Fraction as F

T_WINDOW = 240
FINITE_EXCEPTIONS = [9, 10, 11, 16, 17, 23]   # D-certificate covered


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


def step(p, g, K, alpha=F(1, 4)):
    D = []
    acc = F(0)
    best = None
    for j, pj in enumerate(p):
        acc += pj
        d = acc - F(j + 1, K)
        D.append(d)
        if best is None or d > D[best]:
            best = j
    s = best + 1
    F_s = sum(p[:s])
    M = alpha * F_s
    out = []
    for j, (pj, gj) in enumerate(zip(p, g)):
        base = (1 - alpha) * pj if j < s else pj
        out.append(base + M * gj)
    return out, s, D[best]


def analyze(K: int):
    p0 = bins_p0(K)
    g = bins_g(K)
    m = -(-5 * K // 6)
    b = F(m, K)
    c = 2 * b - b * b
    q = F(6 - K % 6, 3 * K)
    h = F(1, K)
    L = c * q
    states = [p0]
    sels = []
    for _ in range(T_WINDOW):
        p2, s, _ = step(states[-1], g, K)
        sels.append(s)
        states.append(p2)
    # 1. the theorem itself on the trajectory
    U = [K * sum((y - F(1, K)) ** 2 for y in st) for st in states]
    for t in range(T_WINDOW + 1):
        if U[t] < U[2]:
            return (f"THEOREM FAILS: U({t})={U[t]} < U(2)={U[2]}")
    # 2. strengthened decomp, componentwise
    x2 = states[2][m - 1]
    gam = lambda x: (c - x) / (1 - q)
    e = [states[2][i] - gam(x2) * g[i] for i in range(m - 1)]
    for t in range(2, min(60, T_WINDOW)):
        x = states[t][m - 1]
        z = F(3, 4) ** (t - 2)
        for i in range(m - 1):
            if states[t][i] != gam(x) * g[i] + z * e[i]:
                return (f"DECOMP FAILS: t={t} i={i}")
    # 3. DIRECT BACK-SUBSTITUTION (review-10 mandate): the assembled
    #    quadratic Q + 2 z Phi + K E2 z^2 must equal the TRUE U at
    #    t = 2, 3, 4 (both branches occur across the sample)
    G1 = sum(g[:m - 1])
    G2 = sum(gi * gi for gi in g[:m - 1])
    E1 = sum(e)
    GE = sum(g[i] * e[i] for i in range(m - 1))
    E2 = sum(ei * ei for ei in e)
    Sfx = K * sum((p0[i] - F(1, K)) ** 2 for i in range(m, K))
    Phi = lambda x: K * gam(x) * GE - E1
    Qf = lambda x: (K * (gam(x) ** 2 * G2 - 2 * gam(x) * G1 / K
                        + F(m - 1, K * K)) + K * (x - F(1, K)) ** 2 + Sfx)
    for t in (2, 3, 4):
        x = states[t][m - 1]
        z = F(3, 4) ** (t - 2)
        quad = Qf(x) + 2 * z * Phi(x) + K * E2 * z * z
        if quad != U[t]:
            return (f"IDENTITY FAILS: t={t} quad={quad} != U={U[t]} "
                    f"(factor or assembly error)")
    # 4. square-completion sufficient conditions S1-S4 (owner route):
    #    U-U2 = A x^2 + B x + (C-U2) + 2 z (d x + f) + E z^2
    #         = A (x + (B+2dz)/(2A))^2 + H(z);  certify on [0, 3/4].
    #    For the finite exception set these legitimately fail (that is
    #    WHY those K are D-certificate covered); there we check only the
    #    theorem/decomp/identity layers above.
    if K in FINITE_EXCEPTIONS:
        return None
    gA, gB = c / (1 - q), 1 / (1 - q)
    A = K * gB * gB * G2 + K
    B = -2 * K * gA * gB * G2 + 2 * gB * G1 - 2
    Cc = (K * gA * gA * G2 - 2 * gA * G1 + F(m - 1, K) + F(1, K) + Sfx)
    d = -K * gB * GE
    f_ = K * gA * GE - E1
    E = K * E2
    # consistency of the coefficient assembly with Qf/Phi (3 points each)
    for xt in (F(0), F(1, 3), F(5, 7)):
        if A * xt * xt + B * xt + Cc != Qf(xt):
            return "ASSEMBLY FAILS: Q coefficients"
        if d * xt + f_ != Phi(xt):
            return "ASSEMBLY FAILS: Phi coefficients"
    H2 = E - d * d / A
    H1 = 2 * f_ - B * d / A
    Hc = Cc - U[2] - B * B / (4 * A)
    if not A > 0:
        return f"S1 FAILS: A={A}"
    if not H2 > 0:
        return f"S2 FAILS: H2={H2}"
    Hp34 = F(3, 2) * H2 + H1
    if not Hp34 < 0:
        return f"S3 FAILS: H'(3/4)={Hp34}"
    H34 = F(9, 16) * H2 + F(3, 4) * H1 + Hc
    if not H34 >= 0:
        return f"S4 FAILS: H(3/4)={H34}"
    # 5. prefix claims C0
    if not (U[0] >= U[2] and U[1] >= U[2]):
        return "C0 FAILS"
    return None


def main() -> int:
    # sample: every sub-family start (u0 and u0+1 and u0+9), all CYCLE K
    # below 200, the finite difference set, and anchors
    ks = set()
    for r, n_min in ((3, 2), (4, 3), (5, 4)):
        for nu in range(32):
            u0 = 0 if nu >= n_min else 1
            for u in (u0, u0 + 1, u0 + 9):
                n = 32 * u + nu
                ks.add(6 * n + r)
    for K in range(9, 200):
        if K % 6 in (3, 4, 5):
            ks.add(K)
    ks |= {600, 999, 1000, 1001, 2001, 2002, 2003}
    fails = 0
    for K in sorted(ks):
        res = analyze(K)
        if res is not None:
            print(f"[FAIL] K={K}: {res}")
            fails += 1
    print(f"INDEPENDENT E-VERIFY (review-10 repaired): {len(ks)} sampled K "
          f"(window {T_WINDOW}), failures = {fails}")
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
