# -*- coding: utf-8 -*-
"""v04 interval-certifier symbolic kernel (review round-2 I1 rebuild).

Arithmetic PORTED from the v0.3.0 audited implementation
experiments/v03_c_alpha_endpoints/run_c2.py (squarefree / strip_root_at /
count_roots_open / sturm_chain / sign_at; the sturm_selftest there) --
do not fork the semantics.  Added on top:

  pos_closed(a, lo, hi)     a > 0 on the CLOSED [lo, hi]
  nonneg_closed(a, lo, hi)  a >= 0 on [lo, hi] (no interior zeros; only
                            endpoint zeros allowed, verified by sampling)
  pos_leftopen(a, lo, hi)   a > 0 on the LEFT-OPEN (lo, hi]
  nonneg_leftopen(a, lo,hi) a >= 0 on (lo, hi]
  neg_*                     sign-flipped versions

Semantics notes (I1):
  - root counts are DISTINCT real roots in the OPEN (lo, hi);
  - multiple roots: squarefree part first (quotient, remainder asserted
    zero -- the review-I1 bug took the remainder, which is identically
    zero after exact division, collapsing the polynomial);
  - endpoint roots are stripped, so a root exactly at lo/hi does NOT
    count (and never hides a sign change across the open interval);
  - degenerate inputs (zero polynomial, constants) NEVER certify
    positivity -- callers must treat that as OPEN.
"""
from __future__ import annotations

import sys
from fractions import Fraction as F

# deep exact bisection for near-degenerate roots nests ~15 helper frames
# per split level; the split budget below is the real bound
sys.setrecursionlimit(200000)


def ptrim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def peval(a, x):
    r = F(0)
    for c in reversed(a):
        r = r * x + c
    return r


def pdivmod(a, b):
    a = ptrim([F(c) if isinstance(c, int) else c for c in a])
    b = ptrim([F(c) if isinstance(c, int) else c for c in b])
    assert b, "division by zero polynomial"
    db = len(b) - 1
    q = [F(0)] * max(0, len(a) - db)
    while ptrim(a) and len(a) - 1 >= db:
        shift = len(a) - 1 - db
        f = a[-1] / b[-1]
        q[shift] = f
        for i, c in enumerate(b):
            a[shift + i] -= f * c
        a = ptrim(a)
    return ptrim(q), ptrim(a)


def pgcd(a, b):
    while b:
        a, b = b, pdivmod(a, b)[1]
    return a


def pderiv(a):
    return ptrim([a[i] * i for i in range(1, len(a))])


def squarefree(P):
    """Square-free part P / gcd(P, P') (exact; PORT of run_c2)."""
    P = ptrim(list(P))
    if len(P) <= 1:
        return P
    g = pgcd(P, pderiv(P))
    if len(g) <= 1:
        return P
    q, r = pdivmod(P, g)
    assert not r, "squarefree: inexact division"
    return q


def strip_root_at(P, x0):
    """If P(x0) == 0, divide out the factor (x - x0) (PORT)."""
    if peval(P, x0) != 0:
        return P
    q, r = pdivmod(P, [-x0, F(1)])
    assert not r, "strip_root_at: inexact division"
    return q


def sturm_chain(P):
    P = ptrim(list(P))
    if len(P) <= 1:
        return []
    seq = [P, pderiv(P)]
    while True:
        _, r = pdivmod(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-c for c in r])
    return seq


def sign_at(seq, x):
    signs = []
    for p in seq:
        v = peval(p, x)
        if v != 0:
            signs.append(1 if v > 0 else -1)
    changes = 0
    for i in range(len(signs) - 1):
        if signs[i] != signs[i + 1]:
            changes += 1
    return changes


def count_roots_open(P, lo, hi):
    """#distinct roots of P in the OPEN (lo, hi), exact (PORT)."""
    P = ptrim(list(P))
    if not P:
        raise ValueError("count_roots_open: zero polynomial")
    P = squarefree(P)
    if len(P) <= 1:
        return 0
    P = strip_root_at(P, lo)
    P = strip_root_at(P, hi)
    if len(P) <= 1:
        return 0
    seq = sturm_chain(P)
    if not seq:
        return 0
    assert peval(P, lo) != 0 and peval(P, hi) != 0
    return sign_at(seq, lo) - sign_at(seq, hi)


# ---------- root brackets (no width tolerance; exact Sturm counts) ------
class BudgetExceeded(Exception):
    pass


def _isolated_brackets(P, lo, hi, budget=None):
    """Distinct roots of P in the OPEN (lo, hi) as sorted pairwise
    disjoint closed brackets, each containing exactly one distinct root.

    budget: fresh per TOP-LEVEL call (round-3 fix 3.2: the previous
    mutable default [2000] was shared across calls, depleting silently).
    Recursive subcalls share the caller's budget object deliberately.

    Construction (review round-3 fix): recursive midpoint splitting with
    exact counts only -- NO width tolerance.  A midpoint that is exactly
    a root becomes its own bracket (m, m) and both sides are split
    further; therefore any two adjacent brackets are separated by at
    least one split point that is NOT a root (the separator probes below
    rely on this)."""
    budget = [2000] if budget is None else budget
    out = []

    def split(a, b, n, depth):
        if budget[0] <= 0:
            raise BudgetExceeded()
        budget[0] -= 1
        m = (a + b) / 2
        if peval(P, m) == 0:
            out.append((m, m))
            nl = count_roots_open(P, a, m)
            nr = n - 1 - nl
            if nl == 1:
                out.append((a, m))
            elif nl > 1:
                split(a, m, nl, depth + 1)
            if nr == 1:
                out.append((m, b))
            elif nr > 1:
                split(m, b, nr, depth + 1)
            return
        cl = count_roots_open(P, a, m)
        cr = n - cl
        if cl == 0:
            if cr == 1:
                out.append((m, b))
            elif cr > 1:
                split(m, b, cr, depth + 1)
            return
        if cr == 0:
            if cl == 1:
                out.append((a, m))
            elif cl > 1:
                split(a, m, cl, depth + 1)
            return
        if cl == 1:
            out.append((a, m))
        else:
            split(a, m, cl, depth + 1)
        if cr == 1:
            out.append((m, b))
        else:
            split(m, b, cr, depth + 1)

    n = count_roots_open(P, lo, hi)
    if n:
        split(lo, hi, n, 0)
    return sorted(out)


def _probe_left_gap(P, lo, v, budget):
    """Rational point strictly inside (lo, r) where r is the (unique)
    root inside the bracket starting at lo.  Exact bisection on counts;
    terminates because r > lo and dyadics are dense."""
    while budget[0] > 0:
        budget[0] -= 1
        m = (lo + v) / 2
        if peval(P, m) == 0:
            return (lo + m) / 2          # m IS the root; (lo, m) is root-free
        if count_roots_open(P, lo, m) == 0:
            return m                     # m strictly between lo and the root
        v = m
    raise BudgetExceeded()


def _probe_right_gap(P, u, hi, budget):
    """Rational point strictly inside (r, hi) for the root inside the
    bracket ending at hi."""
    while budget[0] > 0:
        budget[0] -= 1
        m = (u + hi) / 2
        if peval(P, m) == 0:
            return (m + hi) / 2
        if count_roots_open(P, m, hi) == 0:
            return m
        u = m
    raise BudgetExceeded()


def _narrow_bracket(P, u, v, budget, side):
    """Halve the bracket [u, v] (exactly one distinct root inside),
    keeping the root; side 'l' pulls the right end in, 'r' the left."""
    if budget[0] <= 0:
        raise BudgetExceeded()
    budget[0] -= 1
    m = (u + v) / 2
    if peval(P, m) == 0:
        return (u, m) if side == "l" else (m, v)
    if count_roots_open(P, u, m) >= 1:
        return (u, m)
    return (m, v)


def _probe_between(P, A, B, budget):
    """Rational point strictly between the two roots of adjacent
    brackets A, B.  Correct regardless of bracket widths: narrows both
    brackets toward their roots; since the roots are distinct and no
    root lies strictly between them, a non-root midpoint separates in
    finitely many dyadic steps."""
    u1, v1 = A
    u2, v2 = B
    while budget[0] > 0:
        budget[0] -= 1
        if v1 < u2:
            mid = (v1 + u2) / 2
            if peval(P, mid) != 0:
                return mid
        u1, v1 = _narrow_bracket(P, u1, v1, budget, "l")
        u2, v2 = _narrow_bracket(P, u2, v2, budget, "r")
    raise BudgetExceeded()


def _gap_probes(P, lo, hi):
    """One exact rational probe strictly inside EVERY connected
    component of (lo, hi) minus the distinct roots of P."""
    brackets = _isolated_brackets(P, lo, hi)
    if not brackets:
        return [(lo + hi) / 2]
    budget = [40000]
    probes = []
    u1, v1 = brackets[0]
    probes.append(_probe_left_gap(P, lo, v1, budget)
                  if u1 == lo else (lo + u1) / 2)
    for A, B in zip(brackets[:-1], brackets[1:]):
        probes.append(_probe_between(P, A, B, budget))
    ul, vl = brackets[-1]
    probes.append(_probe_right_gap(P, ul, hi, budget)
                  if vl == hi else (vl + hi) / 2)
    return probes


# ---------- positivity primitives ----------------------------------------
def pos_closed(a, lo, hi):
    """a > 0 everywhere on [lo, hi].  Degenerate -> False (caller opens)."""
    a = ptrim(list(a))
    if not a:
        return False
    if peval(a, lo) <= 0 or peval(a, hi) <= 0:
        return False
    return count_roots_open(a, lo, hi) == 0


def nonneg_closed(a, lo, hi):
    """a >= 0 on [lo, hi] -- EXACT weak nonnegativity (round-3 rebuild).

    Sign is constant on each component between consecutive distinct
    roots, so it suffices to probe EVERY component with a strictly
    interior rational point (see _gap_probes; no isolation width, no
    fixed tolerance).  Budget exhaustion raises BudgetExceeded -> the
    caller must report UNKNOWN/OPEN."""
    a = ptrim(list(a))
    if not a:
        return True          # zero polynomial is >= 0 (callers decide)
    if peval(a, lo) < 0 or peval(a, hi) < 0:
        return False
    probes = _gap_probes(a, lo, hi)
    return all(peval(a, x) >= 0 for x in probes)


def pos_leftopen(a, lo, hi):
    """a > 0 on the LEFT-OPEN (lo, hi] (value at lo itself irrelevant).
    No roots in (lo, hi) + positive at hi  =>  sign constant '+' on
    (lo, hi] (continuity to the left of hi)."""
    a = ptrim(list(a))
    if not a:
        return False
    if peval(a, hi) <= 0:
        return False
    return count_roots_open(a, lo, hi) == 0


def nonneg_leftopen(a, lo, hi):
    """a >= 0 on (lo, hi] -- exact; the lo-side component IS probed
    (the value at lo itself is irrelevant, but the open gap next to it
    must be nonnegative)."""
    a = ptrim(list(a))
    if not a:
        return True
    if peval(a, hi) < 0:
        return False
    probes = _gap_probes(a, lo, hi)
    return all(peval(a, x) >= 0 for x in probes)


def neg_closed(a, lo, hi):
    return pos_closed([-c for c in a], lo, hi)


def neg_leftopen(a, lo, hi):
    return pos_leftopen([-c for c in a], lo, hi)


# ---------- self-test (I1 boundary cases; review polynomial included) ----
def sturm_selftest():
    ok = True

    def check(name, cond):
        nonlocal ok
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
        ok = ok and cond

    a0 = F(1, 4)
    # review I1 polynomial: (a-1/4)^2 * ((a-1/4)^2 - 1/1600)
    sq = [a0 * a0, -2 * a0, 1]
    sqm = [a0 * a0 - F(1, 1600), -2 * a0, 1]

    def pmul(a, b):
        out = [F(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return ptrim(out)

    P = pmul(sq, sqm)
    lo, hi = F(1, 5), F(3, 10)
    check("I1: P(19/80) < 0 (review value -3/40960000)",
          peval(P, F(19, 80)) == F(-3, 40960000))
    check("I1: count_roots_open = 3 (9/40, 1/4, 11/40)",
          count_roots_open(P, lo, hi) == 3)
    check("I1: pos_closed rejects", not pos_closed(P, lo, hi))
    check("I1: nonneg_closed rejects", not nonneg_closed(P, lo, hi))
    # simple root, sign change
    Q = [-a0, 1]
    check("simple root: count 1 on (0, 1/2)",
          count_roots_open(Q, F(0), F(1, 2)) == 1)
    check("simple root: pos_closed on [0.3, 0.4] accepts",
          pos_closed(Q, F(3, 10), F(2, 5)))
    # double root, no sign change
    R = sq                                     # (a-1/4)^2
    check("double root: count 1 (distinct) on (0,1/2)",
          count_roots_open(R, F(0), F(1, 2)) == 1)
    check("double root: pos_closed on [0, 1/2] rejects (root inside)",
          not pos_closed(R, F(0), F(1, 2)))
    check("double root: nonneg_closed on [0, 1/2] accepts (tangent)",
          nonneg_closed(R, F(0), F(1, 2)))
    # endpoint roots
    check("root at left endpoint not counted in open",
          count_roots_open(Q, a0, F(1, 2)) == 0)
    check("root at right endpoint not counted in open",
          count_roots_open(Q, F(0), a0) == 0)
    check("root at endpoint: pos_closed rejects (value 0 at endpoint)",
          not pos_closed(Q, F(0), a0))
    # zero polynomial / constants
    check("zero polynomial: pos_closed False",
          not pos_closed([], F(0), F(1)))
    check("constant positive accepts", pos_closed([F(1)], F(0), F(1)))
    check("constant negative rejects", not pos_closed([F(-1)], F(0), F(1)))
    try:
        count_roots_open([], F(0), F(1))
        check("zero polynomial raises", False)
    except ValueError:
        check("zero polynomial raises", True)
    # multiple sign changes
    S = pmul([-F(1, 5), 1], pmul([-F(2, 5), 1], [-F(3, 5), 1]))
    check("three sign changes: count 3 on (0,1)",
          count_roots_open(S, F(0), F(1)) == 3)
    check("three sign changes: pos_closed rejects",
          not pos_closed(S, F(0), F(1)))
    # left-open semantics
    check("left-open: (a-1/4) > 0 on (1/4, 1/2]",
          pos_leftopen(Q, a0, F(1, 2)))
    check("left-open: (a-1/4) not > 0 on (0, 1/4]",
          not pos_leftopen(Q, F(0), a0))
    check("left-open: tangent square >= 0 on (1/4, 1/2]",
          nonneg_leftopen(R, a0, F(1, 2)))
    # ---- round-3: near-endpoint / near-double-root exactness -----------
    eps = F(1, 10 ** 40)
    # review counterexample: root AT lo + tiny negative dip inside
    Pn = [a0 * (a0 + eps), -(2 * a0 + eps), F(1)]
    hin = F(2001, 8000)
    check("R3 near-left double root: P(a0+eps/2) < 0",
          peval(Pn, a0 + eps / 2) == -eps * eps / 4)
    check("R3 near-left: nonneg_closed REJECTS",
          not nonneg_closed(Pn, a0, hin))
    check("R3 near-left: nonneg_leftopen REJECTS",
          not nonneg_leftopen(Pn, a0, hin))
    # right-endpoint mirror: dip just left of hi
    Pm = [(hin - eps) * hin, -(2 * hin - eps), F(1)]
    check("R3 near-right mirror: nonneg_closed REJECTS",
          not nonneg_closed(Pm, a0, hin))
    # interior near-double root with negative dip between
    mid = (a0 + hin) / 2
    Pi = [mid * mid - eps * eps / 4, -2 * mid, 1]
    check("R3 interior near-double: nonneg_closed REJECTS",
          not nonneg_closed(Pi, a0, hin))
    check("R3 interior near-double: gap probe finds the dip",
          peval(Pi, mid) == -eps * eps / 4)
    # positive tangential double root still accepted (weak nonneg)
    T = pmul([-(mid), 1], [-(mid), 1])
    check("R3 positive tangential double root: nonneg accepts",
          nonneg_closed(T, a0, hin))
    # strictly positive never accepts zeros
    check("R3 tangential: pos_closed still rejects",
          not pos_closed(T, a0, hin))
    # budget exhaustion raises BudgetExceeded (UNKNOWN path)
    try:
        _isolated_brackets(Pn, a0, hin, budget=[0])
        check("budget exhaustion raises BudgetExceeded", False)
    except BudgetExceeded:
        check("budget exhaustion raises BudgetExceeded", True)
    print("SYMKERNEL SELFTEST:", "ALL PASS" if ok else "FAILURES")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if sturm_selftest() else 1)
