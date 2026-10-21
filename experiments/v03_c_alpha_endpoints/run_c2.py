# -*- coding: utf-8 -*-
"""Research direction C (2026-10-19 handover): general-alpha endpoint
classification via the (K, rho) polynomial template.

Theorem G-ep (K=10, rho=2/5) showed the LN Thm 3(iii) endpoint exception
really occurs at an algebraic alpha: delta_K(alpha) = delta^+(alpha) at the
unique root of a quartic.  Here the template is generalized: for each
(K, rho = p/q) the difference

    diff^+(alpha) = delta_K(alpha) - delta^+(lambda, mu, rho)
    diff^-(alpha) = delta_K(alpha) - delta^-(lambda, mu, rho)

with lambda = 1-alpha, mu = (1-alpha q_K)/(1-alpha), delta_K = 1-(1-alpha)
eta_K (affine in alpha), and delta^+/- the LN Theorem 1 endpoints (S is a
polynomial in lambda, mu with p/q-determined coefficients), is a rational
function of alpha; the numerator is a polynomial with rational
coefficients.  Roots in (0,1) are counted EXACTLY with Sturm chains over
fractions.Fraction, and isolated to rational brackets by exact Sturm
bisection.  Poles (zeros of the denominator) are counted and excluded.

Outputs per (K, rho): N^+ = #{alpha in (0,1): delta_K = delta^+} (endpoint-
exception alphas: at these LN gives rotation number rho but NO periodic
orbit), N^- similarly for the left endpoint (plateau boundary crossings,
no exception), and isolating brackets for every root.

Anchor: (K=10, rho=2/5) must reproduce exactly one root in (0,1) inside
the registered bracket 0.245996 < alpha* < 0.245997 (Theorem G-ep).
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from fractions import Fraction as F
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "v03_p_plateau_allk"))

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(1000000)

from run_a import (RF, pp_add, pp_eval, pp_gcd, pp_mul, ptrim,  # noqa: E402
                   pp_pdivmod, rf_const)
import run_a  # noqa: E402
from sand_m0.word_cyclicity import boundary_params  # noqa: E402

run_a.CANCEL = True   # gcd-cancel rational functions (degrees stay small)

RESULTS = Path(__file__).resolve().parent / "results"

# (K, rho) study set: canonical rho of each K (wp1 table), the K=10
# G-ep case, and neighbours walked in alpha (G sec 5 ladder samples).
TARGETS = [
    (9, 1, 12), (9, 1, 19), (9, 1, 7),
    (10, 7, 18), (10, 9, 22), (10, 2, 5), (10, 1, 2),
    (11, 11, 16),
    (21, 1, 16), (21, 9, 22), (21, 1, 6),
    (65, 11, 16),
    (99, 3, 29),
]


# ---------------- polynomial helpers ----------------
def pp_pdivmod(a, b):
    """(quotient, remainder) of polynomial division, Fractions (inputs
    normalized so bare ints never trigger float division)."""
    a = ptrim([F(c) for c in a])
    b = ptrim([F(c) for c in b])
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


def pp_deriv(a):
    return ptrim([i * c for i, c in enumerate(a)][1:] or [F(0)])


def pp_scale(a, s):
    return [c * s for c in a]


def squarefree(P):
    """Square-free part P / gcd(P, P') (exact)."""
    P = ptrim(list(P))
    if len(P) <= 1:
        return P
    g = pp_gcd(P, pp_deriv(P))
    if len(g) <= 1:
        return P
    q, r = pp_pdivmod(P, g)
    assert not r, "squarefree: inexact division"
    return q


def strip_root_at(P, x0):
    """If P(x0) == 0, divide out the (simple) factor (x - x0)."""
    if pp_eval(P, x0) != 0:
        return P
    q, r = pp_pdivmod(P, [-x0, F(1)])
    assert not r, "strip_root_at: inexact division"
    return q


def sturm_chain(P):
    P = ptrim(list(P))
    if len(P) <= 1:
        return []
    seq = [P, pp_deriv(P)]
    while True:
        _, r = pp_pdivmod(seq[-2], seq[-1])
        if not r:
            break
        seq.append(pp_scale(r, F(-1)))
    return seq


def sign_at(seq, x):
    signs = []
    for p in seq:
        v = pp_eval(p, x)
        if v != 0:
            signs.append(1 if v > 0 else -1)
    changes = 0
    for i in range(len(signs) - 1):
        if signs[i] != signs[i + 1]:
            changes += 1
    return changes


def count_roots_open(P, lo, hi):
    """#roots of P in the OPEN interval (lo, hi), exact.  Review-7 fix:
    square-free decomposition first (multiplicity collapses, so multiple
    roots are counted once), then roots sitting exactly at the endpoints
    are divided out (they are excluded from the open interval), then plain
    Sturm with P(lo), P(hi) != 0.  The zero polynomial is an error (the
    difference would be identically zero -- infinitely many solutions)."""
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
    assert pp_eval(P, lo) != 0 and pp_eval(P, hi) != 0
    return sign_at(seq, lo) - sign_at(seq, hi)   # roots in (lo, hi)


def bisect_root(P, lo, hi, depth=200):
    """Isolate ONE root known to be in (lo, hi) (with P(lo), P(hi) != 0 and
    exactly one root inside) to a narrow rational bracket."""
    for _ in range(depth):
        mid = (lo + hi) / 2
        if pp_eval(P, mid) == 0:
            return (mid, mid)
        # sturm decides which half keeps the root
        if count_roots_open(P, lo, mid) == 1:
            hi = mid
        else:
            lo = mid
        if hi - lo < F(1, 10 ** 30):
            break
    return (lo, hi)


def isolate_all(P, lo=F(0), hi=F(1)):
    """All roots of P in (lo, hi) with exact rational brackets."""
    total = count_roots_open(P, lo, hi)
    if total == 0:
        return []
    out = []
    frontier = [(lo, hi, total)]
    while frontier:
        a, b, n = frontier.pop()
        if n == 1:
            out.append(bisect_root(P, a, b))
            continue
        m = (a + b) / 2
        if pp_eval(P, m) == 0:
            # root exactly at m: isolate to a point and split around it
            out.append((m, m))
            nl = count_roots_open(P, a, m)
            nr = n - 1 - nl
            if nl:
                frontier.append((a, m, nl))
            if nr:
                frontier.append((m, b, nr))
            continue
        nl = count_roots_open(P, a, m)
        nr = count_roots_open(P, m, b)
        assert nl + nr == n, "sturm split mismatch"
        if nl:
            frontier.append((a, m, nl))
        if nr:
            frontier.append((m, b, nr))
    return sorted(out)


# ---------------- symbolic (K, rho) differences in alpha ----------------
def diffs_rf(K, p, q):
    """(diff_plus, diff_minus) as rational functions of alpha."""
    P = boundary_params(K)
    qK, c, h = P["q"], P["c"], P["h"]
    eta = qK * (c - h) / ((1 - qK) * h)
    a = RF([F(0), F(1)], [F(1)])           # alpha
    lam = rf_const(1) - a
    mu = (rf_const(1) - a * qK) / lam
    deltaK = rf_const(1) - lam * eta
    # S = sum_{k=1}^{q-2} d_k lambda^k mu^{floor(kp/q)}
    S = rf_const(0)
    for k in range(1, q - 1):
        d = ((k + 1) * p) // q - (k * p) // q
        S = S + (lam ** k) * (mu ** ((k * p) // q)) * d
    one = rf_const(1)
    nump = (one - lam) * (one + mu * S + (lam ** (q - 1)) * (mu ** p) *
                          (one - lam))
    denp = one + (mu - one) * S + (lam ** (q - 1)) * (mu ** (p - 1)) * (
        mu - lam * mu - one)
    numm = (one - lam) * (one + mu * S)
    denm = one + (mu - one) * S - (lam ** q) * (mu ** (p - 1))
    return deltaK - nump / denp, deltaK - numm / denm


def classify(K, p, q):
    dp, dm = diffs_rf(K, p, q)
    res = {}
    for tag, diff in (("plus", dp), ("minus", dm)):
        N, D = diff.num, diff.den
        n_roots = count_roots_open(N, F(0), F(1))
        n_poles = count_roots_open(D, F(0), F(1))
        brackets = []
        if n_roots:
            for (lo, hi) in isolate_all(N, F(0), F(1)):
                # a common root of N and D would be a removable artifact,
                # not a genuine crossing; probe both at the bracket point
                if lo == hi:
                    genuine = pp_eval(N, lo) == 0 and pp_eval(D, lo) != 0
                else:
                    genuine = True
                brackets.append({"bracket": [str(lo), str(hi)],
                                 "genuine": genuine})
        res[tag] = {"num_degree": len(N) - 1, "den_degree": len(D) - 1,
                    "roots_open_0_1": n_roots, "poles_open_0_1": n_poles,
                    "brackets": brackets}
    return res


def sturm_selftest():
    """Exact counterexamples from owner review 7 (multiple roots and
    endpoint roots): the old counter returned 0 for x^2(x-1/2) on (0,1)."""
    cases = [
        ([0, 0, F(-1, 2), 1], 1),          # x^2 (x - 1/2): root 1/2 only
        ([F(1, 4), F(-1), 1], 1),          # (x - 1/2)^2: double root counted once
        ([0, 0, 1], 0),                    # x^2: root only AT endpoint 0
        ([1, -2, 1], 0),                   # (x-1)^2: root only AT endpoint 1
        ([F(-1, 2), F(1, 2), -2, 1], 1),   # (x-1)^2 (x-1/2)... see below
        ([F(5)], 0),                       # constant
        ([0, 1], 0),                       # x: root at endpoint only
    ]
    # (x-1)^2 (x-1/2) = (x^2-2x+1)(x-1/2) = x^3 - 2.5x^2 + 2x - 0.5
    cases[4] = ([F(-1, 2), 2, F(-5, 2), 1], 1)
    for P, want in cases:
        got = count_roots_open(P, F(0), F(1))
        assert got == want, (P, got, want)
    try:
        count_roots_open([], F(0), F(1))
        raise AssertionError("zero polynomial must raise")
    except ValueError:
        pass
    return True


def verify_table(data):
    """COMPLETE acceptance for a C endpoint-classification table (owner
    review 9): the single source of truth -- the generator self-check and
    every positive/negative test call THIS function.  Checks, per side
    (plus/minus) of every target:
      * the target set equals TARGETS exactly, with no duplicates;
      * root count, pole count, numerator AND denominator degrees match a
        fresh recompute (diffs_rf rebuilt from model constants);
      * bracket count equals the root count; every bracket lies in (0,1);
      * each bracket isolates EXACTLY ONE numerator root (Sturm), with no
        pole inside, and genuine=True;
      * the brackets JOINTLY isolate the target roots exactly once (no
        gaps, no double-isolation of the same root);
      * zero-width point brackets are accepted only with an exact
        numerator zero at the point and a nonzero denominator there;
      * the source hash binds to the running run_c2.py.
    Returns a failure list (empty = accepted)."""
    import hashlib
    fails = []
    sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if data.get("source_sha256") != sha:
        fails.append("source_sha256 mismatch")
    got = [(t["K"], int(t["rho"].split("/")[0]), int(t["rho"].split("/")[1]))
           for t in data.get("targets", [])]
    if sorted(got) != sorted(TARGETS) or len(set(got)) != len(got):
        fails.append(f"target set mismatch or duplicates: {sorted(got)}")
        return fails
    for t in data["targets"]:
        p, q = (int(x) for x in t["rho"].split("/"))
        K = t["K"]
        dp, dm = diffs_rf(K, p, q)
        for side, diff in (("plus", dp), ("minus", dm)):
            num, den = diff.num, diff.den
            rec = t[side]
            for field, want in (
                    ("roots_open_0_1", count_roots_open(num, F(0), F(1))),
                    ("poles_open_0_1", count_roots_open(den, F(0), F(1))),
                    ("num_degree", len(num) - 1),
                    ("den_degree", len(den) - 1)):
                if rec.get(field) != want:
                    fails.append(f"K={K} rho={p}/{q} {side} {field}: "
                                 f"stored {rec.get(field)} != {want}")
            n_roots = count_roots_open(num, F(0), F(1))
            brs = sorted(rec.get("brackets", []),
                         key=lambda b: F(b["bracket"][0]))
            if len(brs) != n_roots:
                fails.append(f"K={K} rho={p}/{q} {side}: {len(brs)} "
                             f"brackets != {n_roots} roots")
                continue
            isolated = 0
            for br in brs:
                lo, hi = F(br["bracket"][0]), F(br["bracket"][1])
                if not (0 < lo <= hi < 1):
                    fails.append(f"K={K} {side}: bracket outside (0,1]: "
                                 f"[{lo}, {hi}]")
                    continue
                if not br.get("genuine", False):
                    fails.append(f"K={K} {side}: bracket not genuine")
                if lo == hi:
                    # zero-width convention: exact numerator root at the
                    # point, denominator nonzero there
                    if pp_eval(num, lo) != 0:
                        fails.append(f"K={K} {side}: zero-width bracket "
                                     f"not an exact root")
                    elif pp_eval(den, lo) == 0:
                        fails.append(f"K={K} {side}: zero-width bracket "
                                     f"is also a pole")
                    isolated += 1
                else:
                    c = count_roots_open(num, lo, hi)
                    if c != 1:
                        fails.append(f"K={K} {side}: bracket "
                                     f"[{lo},{hi}] isolates {c} roots")
                    isolated += c
                    if count_roots_open(den, lo, hi) != 0:
                        fails.append(f"K={K} {side}: pole inside bracket")
            if isolated != n_roots:
                fails.append(f"K={K} rho={p}/{q} {side}: brackets jointly "
                             f"isolate {isolated} != {n_roots} roots "
                             f"(gap or double-isolation)")
    return fails


def main():
    t0 = time.time()
    print("== C: general-alpha endpoint classification (Sturm, exact; "
          "review-7 fixed) ==")
    assert sturm_selftest()
    print("[C] Sturm counterexample self-test PASS "
          "(x^2(x-1/2) -> 1, endpoint/multiple roots handled)")
    out = {"targets": []}
    # anchor: Theorem G-ep reproduction
    anchor = classify(10, 2, 5)
    ap = anchor["plus"]
    assert ap["roots_open_0_1"] == 1, ap
    lo, hi = F(ap["brackets"][0]["bracket"][0]), \
        F(ap["brackets"][0]["bracket"][1])
    assert F(245996, 10 ** 6) < hi and lo < F(245997, 10 ** 6), (lo, hi)
    print(f"[C] anchor K=10 rho=2/5: N^+ = 1, alpha* in "
          f"({float(lo):.9f}, {float(hi):.9f})  [G-ep bracket 0.245996/"
          f"0.245997 reproduced]")
    for (K, p, q) in TARGETS:
        res = classify(K, p, q)
        res["K"], res["rho"] = K, f"{p}/{q}"
        out["targets"].append(res)
        print(f"[C] K={K:3d} rho={p}/{q:2d}: "
              f"N^+={res['plus']['roots_open_0_1']:2d} "
              f"(poles {res['plus']['poles_open_0_1']}) "
              f"N^-={res['minus']['roots_open_0_1']:2d} "
              f"(poles {res['minus']['poles_open_0_1']}) "
              f"deg num^+ {res['plus']['num_degree']}")
    out["git"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                                capture_output=True, text=True).stdout.strip()
    import hashlib
    out["source_sha256"] = hashlib.sha256(
        Path(__file__).read_bytes()).hexdigest()
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "c_alpha_endpoints.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")
    fails = verify_table(out)
    if fails:
        raise SystemExit(f"[C] ABORT: table self-verification failed: "
                         f"{fails[:5]}")
    print(f"[C] PASS ({time.time() - t0:.1f}s); table self-verified "
          f"via verify_table")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
