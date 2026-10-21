# -*- coding: utf-8 -*-
"""Research direction A (2026-10-19 handover): all-K analytic plateau
classification for the canonical alpha=1/4 dictionary.

WP1/F3 certified plateau membership lattice point by lattice point
(301/301).  Here the three INFINITE CYCLE zones are upgraded to
polynomial sign certificates in n (K = 6n + r):

    r=3, rho=1/9,   n0=32   (K >= 195, word BA^8)
    r=4, rho=2/5,   n0=2    (K >= 16,  word BABAA)
    r=5, rho=11/16, n0=10   (K >= 65,  word ABBABBABBBABBABB)

For each zone, with exact Fraction polynomial arithmetic:

    LEFT  : delta_K(n) - delta^-(n) >= 0     for all n >= n0  (left-closed)
    RIGHT : delta^+(n) - delta_K(n) > 0      for all n >= n0  (right-strict;
            this is the canonical LN Thm 3(iii) endpoint exclusion, G3 tail)

Each comparison sign(num)*sign(den) uses two certificates (numerator sign,
denominator positivity).  Certificate format: hybrid = explicit exact
finite values on [n0, N) + shift certificate at N (P(N+m) expanded in m,
all coefficients >= 0, and constant term > 0 for strict claims).

The runner GENERATES the certificates, writes them to JSON, then
RE-VERIFIES from the JSON (fresh polynomial rebuild), and cross-checks
the symbolic dictionary against boundary_params at anchor K.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from fractions import Fraction as F
from pathlib import Path

# uncanceled rational-function products produce very large integer
# coefficients; raise CPython's int->str digit limit for certificate I/O
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(1000000)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from sand_m0.word_cyclicity import boundary_params  # noqa: E402

LAM = F(3, 4)
RESULTS = Path(__file__).resolve().parent / "results"
MAX_SHIFT_SCAN = 600

ZONES = [
    {"r": 3, "p": 1, "q": 9, "n0": 32},     # K >= 195, BA^8
    {"r": 4, "p": 2, "q": 5, "n0": 2},      # K >= 16,  BABAA
    {"r": 5, "p": 11, "q": 16, "n0": 10},   # K >= 65,  ABBABBABBBABBABB
]


# ---------------- exact polynomial arithmetic (Fractions) ----------------
def ptrim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def pp_add(a, b):
    out = [F(0)] * max(len(a), len(b))
    for i, c in enumerate(a):
        out[i] += c
    for i, c in enumerate(b):
        out[i] += c
    return ptrim(out)


def pp_neg(a):
    return [-c for c in a]


def pp_mul(a, b):
    if not a or not b:
        return []
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            out[i + j] += x * y
    return ptrim(out)


def pp_pow(a, e):
    out = [F(1)]
    for _ in range(e):
        out = pp_mul(out, a)
    return out


def pp_eval(a, x):
    tot = F(0)
    for c in reversed(a):
        tot = tot * x + c
    return tot


def pp_shift(a, n0):
    """Coefficients of P(n0 + m) in m, exact (repeated synthetic division
    by (x - n0); successive remainders are P^{(j)}(n0)/j!).  Storage is
    low-to-high, so the recurrence runs from the top coefficient down."""
    coeffs = list(a)
    shifted = []
    while coeffs:
        b = [F(0)] * len(coeffs)
        for i in range(len(coeffs) - 1, -1, -1):
            b[i] = coeffs[i] + (n0 * b[i + 1] if i + 1 < len(coeffs) else 0)
        shifted.append(b[0])          # remainder = next Taylor coefficient
        coeffs = b[1:]                # quotient coefficients (low to high)
    return ptrim(shifted)


def pp_mod(a, b):
    """Remainder of polynomial division (Fractions; inputs normalized --
    bare-int coefficients must not leak float division)."""
    a = ptrim([F(c) for c in a])
    b = ptrim([F(c) for c in b])
    assert b, "poly mod by zero"
    db = len(b) - 1
    while a and len(a) - 1 >= db:
        shift = len(a) - 1 - db
        f = a[-1] / b[-1]
        for i, c in enumerate(b):
            a[shift + i] -= f * c
        a = ptrim(a)
    return a


def pp_gcd(a, b):
    """Monic GCD of two polynomials (Euclidean, exact)."""
    a, b = ptrim(list(a)), ptrim(list(b))
    while b:
        a, b = b, pp_mod(a, b)
    if a:
        lead = a[-1]
        a = [c / lead for c in a]
    return a


def pp_divexact(a, g):
    """Exact division of a by g (g divides a); returns quotient."""
    q, r = pp_pdivmod(a, g)
    assert not r, "pp_divexact: remainder not zero"
    return q


def pp_pdivmod(a, b):
    """(quotient, remainder) of polynomial division, Fractions."""
    a = ptrim(list(a))
    b = ptrim(list(b))
    assert b, "division by zero polynomial"
    db = len(b) - 1
    q = [F(0)] * max(0, len(a) - db)
    while a and len(a) - 1 >= db:
        shift = len(a) - 1 - db
        f = a[-1] / b[-1]
        q[shift] = f
        for i, c in enumerate(b):
            a[shift + i] -= f * c
        a = ptrim(a)
    return ptrim(q), ptrim(a)


CANCEL = False   # run_b2 sets True: cancel gcd(num,den) + sign-normalize.
# With CANCEL False the representation is exactly the one used to generate
# the committed a_plateau_cert.json (reproducibility contract).


class RF:
    """Rational function in n: numerator/denominator coefficient lists."""

    def __init__(self, num, den):
        num, den = ptrim(list(num)), ptrim(list(den))
        assert den, "zero denominator"
        if CANCEL:
            g = pp_gcd(num, den)
            if g and len(g) > 1:
                num = pp_divexact(num, g)
                den = pp_divexact(den, g)
            if den and den[-1] < 0:   # canonical positive leading coeff
                num = [-c for c in num]
                den = [-c for c in den]
        self.num = num
        self.den = den

    def __add__(self, o):
        o = o if isinstance(o, RF) else RF([F(o)], [F(1)])
        return RF(pp_add(pp_mul(self.num, o.den), pp_mul(o.num, self.den)),
                  pp_mul(self.den, o.den))

    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, RF) else RF([F(o)], [F(1)])
        return self + RF(pp_neg(o.num), o.den)

    def __rsub__(self, o):
        return (o if isinstance(o, RF) else RF([F(o)], [F(1)])) - self

    def __mul__(self, o):
        o = o if isinstance(o, RF) else RF([F(o)], [F(1)])
        return RF(pp_mul(self.num, o.num), pp_mul(self.den, o.den))

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = o if isinstance(o, RF) else RF([F(o)], [F(1)])
        return RF(pp_mul(self.num, o.den), pp_mul(self.den, o.num))

    def __rtruediv__(self, o):
        o = o if isinstance(o, RF) else RF([F(o)], [F(1)])
        return RF(pp_mul(o.num, self.den), pp_mul(o.den, self.num))

    def pow(self, e):
        return RF(pp_pow(self.num, e), pp_pow(self.den, e))

    def __pow__(self, e):
        assert isinstance(e, int) and e >= 0
        return self.pow(e)


def rf_const(c):
    return RF([F(c)], [F(1)])


def rf_n():
    return RF([F(0), F(1)], [F(1)])


# ---------------- symbolic dictionary along K = 6n + r ----------------
def dictionary_rf(r):
    """Exact rational-function dictionary in n for K = 6n + r, r in
    {3,4,5}: h=1/K, b=m/K with m=5n+r, c=2b-b^2, q=(6-r)/(3K),
    beta=1-q/4, mu=4 beta/3, eta=q(c-h)/((1-q)h), delta_K=1-(3/4)eta."""
    n = rf_n()
    K = n * 6 + r
    m = n * 5 + r                       # ceil(5K/6) on these residues
    b = m / K
    h = rf_const(1) / K
    q = rf_const(6 - r) / (K * 3)
    c = b * 2 - b * b
    beta = rf_const(1) - q * F(1, 4)
    mu = beta * F(4, 3)
    eta = q * (c - h) / (rf_const(1) - q) * K   # division by h = 1/K
    delta = rf_const(1) - eta * F(3, 4)
    return {"K": K, "m": m, "b": b, "h": h, "q": q, "c": c,
            "beta": beta, "mu": mu, "eta": eta, "delta": delta}


def S_sum_rf(mu, p, qd):
    """S = sum_{k=1}^{qd-2} d_k lambda^k mu^{floor(kp/qd)} with
    d_k = floor((k+1)p/qd) - floor(kp/qd); exact constants, RF powers."""
    if qd == 2:
        return rf_const(0)
    tot = rf_const(0)
    for k in range(1, qd - 1):
        d = ((k + 1) * p) // qd - (k * p) // qd
        tot = tot + mu.pow((k * p) // qd) * (d * LAM ** k)
    return tot


def comparisons_rf(r, p, qd):
    """(left, right) rational functions:
    left = delta_K - delta^-  (claim: >= 0),
    right = delta^+ - delta_K (claim: >  0)."""
    D = dictionary_rf(r)
    mu = D["mu"]
    S = S_sum_rf(mu, p, qd)
    one = rf_const(1)
    lamq1 = LAM ** (qd - 1)
    num_plus = (one - LAM) * (one + mu * S + mu.pow(p) * (lamq1 * (1 - LAM)))
    den_plus = one + (mu - one) * S + mu.pow(p - 1) * (
        lamq1 * (mu * (1 - LAM) - one))
    num_minus = (one - LAM) * (one + mu * S)
    den_minus = one + (mu - one) * S - mu.pow(p - 1) * LAM ** qd
    left = D["delta"] - num_minus / den_minus
    right = num_plus / den_plus - D["delta"]
    return D, left, right


# ---------------- certificates ----------------
def shift_ok(P, N, strict):
    co = pp_shift(P, N)
    if not co:
        return None
    if any(c < 0 for c in co):
        return None
    if strict and not co[0] > 0:
        return None
    return [str(c) for c in co]


def hybrid_certificate(P, n0, strict):
    """First N >= n0 where the shift certificate holds; explicit exact
    finite values cover [n0, N).  The stored 'covers_from' binds the
    claimed start: verification must check the finite prefix is the FULL
    consecutive range [covers_from, N)."""
    for N in range(n0, n0 + MAX_SHIFT_SCAN):
        co = shift_ok(P, N, strict)
        if co is not None:
            finite = [{"n": t, "val": str(pp_eval(P, t))}
                      for t in range(n0, N)]
            return {"strict": strict, "N": N, "covers_from": n0,
                    "shift_coeffs": co, "finite": finite}
    return None


def verify_certificate(P, cert, expected_start=None, expected_strict=None):
    """Re-verify a stored certificate against polynomial P (owner review 8
    hardening).  expected_start / expected_strict are the THEOREM's
    requirements supplied by the CALLER (zone n0 / claim direction), never
    taken from the certificate alone -- a self-consistent certificate
    whose start was silently shifted later (start and N moved together,
    coefficients recomputed, prefix emptied) must be rejected.  Checks:
    certificate's covers_from/strict AGREE with the expectation, full
    consecutive finite coverage from the required start to N, exact
    recomputed values with the required direction, shift coefficients
    match and non-negativity."""
    fails = []
    if expected_start is not None and cert.get("covers_from") != expected_start:
        fails.append(f"certificate covers_from {cert.get('covers_from')} != "
                     f"required start {expected_start}")
        return fails
    if expected_strict is not None and cert.get("strict") != expected_strict:
        fails.append(f"certificate strict={cert.get('strict')} != required "
                     f"strict={expected_strict}")
        return fails
    strict = expected_strict if expected_strict is not None else cert["strict"]
    covers_from = (expected_start if expected_start is not None
                   else cert.get("covers_from"))
    if covers_from is None:
        fails.append("certificate does not bind a claimed start "
                     "(covers_from missing)")
        return fails
    N = cert["N"]
    if not (covers_from <= N):
        fails.append("covers_from > N")
        return fails
    # finite prefix must be exactly the consecutive range [covers_from, N)
    expect_range = list(range(covers_from, N))
    got_range = [rec["n"] for rec in cert["finite"]]
    if got_range != expect_range:
        fails.append(f"finite prefix not the full range [{covers_from},"
                     f" {N}) -- got {got_range[:5]}...")
        return fails
    for rec in cert["finite"]:
        v = pp_eval(P, rec["n"])
        if str(v) != rec["val"]:
            fails.append(f"finite value mismatch at n={rec['n']}")
        ok = v > 0 if strict else v >= 0
        if not ok:
            fails.append(f"finite point n={rec['n']} violates claim")
    co = pp_shift(P, N)
    if [str(c) for c in co] != cert["shift_coeffs"]:
        fails.append("shift coefficients mismatch")
    if any(c < 0 for c in co):
        fails.append("negative coefficient in shift")
    if strict and not co[0] > 0:
        fails.append("strict claim but constant term not positive")
    return fails


def _source_sha256():
    import hashlib
    return hashlib.sha256(
        Path(__file__).read_bytes()).hexdigest()


def verify_from_json():
    """Load the committed certificate JSON and re-verify (owner review 7
    hardening): the zone set is EXACTLY the expected three zones (no
    missing/extra/duplicate), every stored certificate re-verifies against
    freshly rebuilt polynomials with full finite-prefix coverage from the
    claimed start, the generator source hash matches the running file, and
    the wp1_rows.csv cross-check holds.  Returns failure list."""
    import csv
    path = RESULTS / "a_plateau_cert.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    fails = []
    if data.get("source_sha256") != _source_sha256():
        fails.append("source_sha256 mismatch: certificate file generated "
                     "by a different version of run_a.py")
    expected = [(z["r"], z["p"], z["q"], z["n0"]) for z in ZONES]
    got = []
    for entry in data["zones"]:
        p = int(entry["rho"].split("/")[0])
        qd = int(entry["rho"].split("/")[1])
        got.append((entry["r"], p, qd, entry["n0"]))
    if sorted(got) != sorted(expected):
        fails.append(f"zone set mismatch: expected {expected}, got {got}")
    wp1 = {int(row["K"]): (int(row["rho_p"]), int(row["rho_q"]))
           for row in csv.DictReader(
               (REPO / "experiments" / "v03_wp1_plateau" / "results" /
                "wp1_rows.csv").open(encoding="utf-8"))}
    for entry in data["zones"]:
        r, p, qd, n0 = entry["r"], None, None, entry["n0"]
        p = int(entry["rho"].split("/")[0])
        qd = int(entry["rho"].split("/")[1])
        if (r, p, qd, n0) not in expected:
            continue   # already reported above
        D, left, right = comparisons_rf(r, p, qd)
        claims = {"left_num": (left.num, False),
                  "left_den": (left.den, True),
                  "right_num": (right.num, True),
                  "right_den": (right.den, True)}
        stored_names = set(entry["certificates"])
        if stored_names != set(claims):
            fails.append(f"r={r}: claim set mismatch {stored_names} != "
                         f"{set(claims)}")
            continue
        for name, (P, strict) in claims.items():
            fails += [f"r={r} {name}: {m}" for m in
                      verify_certificate(P, entry["certificates"][name],
                                         expected_start=n0,
                                         expected_strict=strict)]
        # committed wp1 plateau table must agree on sampled zone points
        # (census lattice 9..603 + anchors 999/1000/1001)
        for K in (6 * n0 + r, 594 + r, 996 + r):
            if wp1.get(K) != (p, qd):
                fails.append(f"r={r}: wp1_rows rho at K={K} is {wp1.get(K)}, "
                             f"expected {(p, qd)}")
    return fails


def main():
    t0 = time.time()
    print("== A: all-K plateau polynomial certificates (canonical a=1/4) ==")
    out = {"zones": [], "config": {"experiment": "v03_p_plateau_allk",
                                   "lambda": "3/4"}}
    for z in ZONES:
        r, p, qd, n0 = z["r"], z["p"], z["q"], z["n0"]
        D, left, right = comparisons_rf(r, p, qd)
        claims = {
            "left_num": (left.num, False),
            "left_den": (left.den, True),
            "right_num": (right.num, True),
            "right_den": (right.den, True),
        }
        certs = {}
        for name, (P, strict) in claims.items():
            c = hybrid_certificate(P, n0, strict)
            if c is None:
                raise SystemExit(
                    f"[A] ABORT r={r} {name}: no shift certificate within "
                    f"+{MAX_SHIFT_SCAN} (degree {len(P)-1})")
            certs[name] = c
        # independent re-verification from scratch
        verification = {name: verify_certificate(P, certs[name])
                        for name, (P, _) in claims.items()}
        if any(verification.values()):
            raise SystemExit(f"[A] ABORT r={r}: verification failed: "
                             f"{verification}")
        # symbolic-vs-numeric anchors (exact): dictionary AND both
        # plateau comparisons against the run_wp1 formulas
        anchors = []
        for K in (6 * n0 + r, 6 * n0 + 6 + r, 240 + r, 600 + r):
            nn = (K - r) // 6
            P = boundary_params(K)
            mu_sym = pp_eval(D["mu"].num, nn) / pp_eval(D["mu"].den, nn)
            d_sym = pp_eval(D["delta"].num, nn) / pp_eval(D["delta"].den, nn)
            mu_ex = 4 * P["beta"] / 3
            eta_ex = P["q"] * (P["c"] - P["h"]) / ((1 - P["q"]) * P["h"])
            d_ex = 1 - F(3, 4) * eta_ex
            S_ex = F(0)
            for k in range(1, qd - 1):
                dk = ((k + 1) * p) // qd - (k * p) // qd
                S_ex += dk * LAM ** k * mu_ex ** ((k * p) // qd)
            dlo = ((1 - LAM) * (1 + mu_ex * S_ex)) / (
                1 + (mu_ex - 1) * S_ex - LAM ** qd * mu_ex ** (p - 1))
            dhi = ((1 - LAM) * (1 + mu_ex * S_ex + LAM ** (qd - 1) *
                                mu_ex ** p * (1 - LAM))) / (
                1 + (mu_ex - 1) * S_ex + LAM ** (qd - 1) * mu_ex ** (p - 1) *
                (mu_ex - LAM * mu_ex - 1))
            lv = pp_eval(left.num, nn) / pp_eval(left.den, nn)
            rv = pp_eval(right.num, nn) / pp_eval(right.den, nn)
            if (mu_sym != mu_ex or d_sym != d_ex or lv != d_ex - dlo
                    or rv != dhi - d_ex):
                raise SystemExit(f"[A] ABORT: symbolic/numeric mismatch K={K}")
            anchors.append({"K": K, "ok": True})
        # plateau membership anchor: numeric at n0 and one beyond
        entry = {
            "r": r, "rho": f"{p}/{qd}", "n0": n0,
            "degrees": {k: len(P) - 1 for k, (P, _) in claims.items()},
            "N_shift": {k: c["N"] for k, c in certs.items()},
            "finite_counts": {k: len(c["finite"]) for k, c in certs.items()},
            "certificates": {k: c for k, c in certs.items()},
            "verification_all_empty": True,
            "anchors": anchors,
        }
        out["zones"].append(entry)
        print(f"[A] r={r} rho={p}/{qd} n0={n0}: degrees {entry['degrees']}"
              f" shift-N {entry['N_shift']} finite {entry['finite_counts']}"
              f" -- VERIFIED")
    out["git"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                                capture_output=True, text=True).stdout.strip()
    out["source_sha256"] = _source_sha256()
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "a_plateau_cert.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")
    fails = verify_from_json()
    if fails:
        raise SystemExit(f"[A] ABORT: JSON re-verification failed: {fails}")
    print(f"[A] PASS ({time.time() - t0:.1f}s); certificates -> "
          f"{RESULTS / 'a_plateau_cert.json'} (+ JSON re-verification, "
          f"wp1_rows cross-check)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
