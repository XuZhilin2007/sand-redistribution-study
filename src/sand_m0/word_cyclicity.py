"""v0.3 word-cyclicity theory: exact eventual-periodicity machinery for the
canonical XO boundary scalar map (alpha = 1/4), building on S3 T1-T5.

Objects
-------
Scalar map (T4, exact): with h = 1/K, c = 2b - b^2, q = (6 - r)/(3K),
L = c q and boundary coordinate x = p_m,

    A(x) = (3/4) x + L/4   (x > h),      B(x) = x + (c - x) q/4  (x <= h).

Phase reduction (Lemma B here): on the A-landing phase band
xi = (x - L)/(h - L) in (3/4, 1] (r = 3, 4: B-runs have length exactly 1)
or on the B-crossing band (r = 5), the map is a piecewise-affine contraction
with explicitly known cells.  For r in {3, 4}

    xi_B = beta xi + sigma,   beta = 1 - q/4,   sigma = q (c - L) / (4 (h - L)),

followed by an A-run of length n = min{k : (3/4)^k xi_B <= 1}.  A pure window
sigma in ((4/3)^(k-1) - (3/4) beta, (4/3)^k - beta] makes n == k on the whole
band (global island: every orbit locks to the word A^k B).

Certificates: for a word w over {A, B} with affine composition
Phi_w(x) = lam x + mu and fixed point x* = mu / (1 - lam), an interval
[x* - delta, x* + delta] whose one-period orbit chain stays strictly on the
branch sides of h and closes up proves that every point of the interval has
itinerary w^infty.

All arithmetic is exact Fraction; no floating point is used in any check.
"""
from __future__ import annotations

import math
from fractions import Fraction

ALPHA = Fraction(1, 4)


def boundary_params(K: int) -> dict:
    """Exact boundary parameters for integer K >= 6 (S3 notation).

    m = ceil(5K/6) is computed in pure integer arithmetic (no float):
    K up to any magnitude stays exact."""
    if K < 6:
        raise ValueError("K >= 6 required")
    m = (5 * K + 5) // 6  # == ceil(5K/6)
    assert 6 * m - 5 <= 5 * K <= 6 * m
    h = Fraction(1, K)
    b = Fraction(m, K)
    c = 2 * b - b * b
    q = Fraction(6 - (K % 6), 3 * K)
    L = c * q
    return {"K": K, "r": K % 6, "m": m, "h": h, "b": b, "c": c,
            "q": q, "L": L, "a": b * (1 - b), "beta": 1 - q * ALPHA,
            "sigma": q * (c - L) / (4 * (h - L)) if h != L else None}


def scalar_step(x: Fraction, P: dict):
    """One exact scalar step; returns (x', branch)."""
    if x > P["h"]:
        return Fraction(3, 4) * x + P["L"] * ALPHA, "A"
    return x + (P["c"] - x) * P["q"] * ALPHA, "B"


def forced_step(x: Fraction, P: dict, branch: str) -> Fraction:
    if branch == "A":
        return Fraction(3, 4) * x + P["L"] * ALPHA
    return x + (P["c"] - x) * P["q"] * ALPHA


def scalar_orbit(x0: Fraction, P: dict, T: int):
    xs, it = [x0], []
    x = x0
    for _ in range(T):
        x, br = scalar_step(x, P)
        xs.append(x)
        it.append(br)
    return xs, it


def runs_of(it):
    out = []
    for ch in it:
        if out and out[-1][0] == ch:
            out[-1][1] += 1
        else:
            out.append([ch, 1])
    return out


def detect_eventual_word(it, rep: int = 20, t0_cap: int = 80,
                         per_cap: int = 80):
    """Strong detection: run sequence (last, possibly truncated run dropped)
    is eventually periodic with >= rep repeats. Returns (i0, n, word)."""
    items = [tuple(v) for v in runs_of(it)][:-1]
    t0_cap = min(t0_cap, max(0, len(items) - rep * 2))
    per_cap = min(per_cap, len(items) // rep)
    for t0 in range(t0_cap + 1):
        tail = items[t0:]
        for per in range(1, per_cap + 1):
            seg = tail[:per]
            if all(tail[i] == seg[i % per] for i in range(len(tail))):
                i0 = sum(ln for _, ln in items[:t0])
                n = sum(ln for _, ln in seg)
                return i0, n, "".join(it[i0:i0 + n])
    return None, None, None


def compose_word(w: str, P: dict):
    """Affine composition of the word's branches: x -> lam x + mu."""
    lam, mu = Fraction(1), Fraction(0)
    for br in w:
        if br == "A":
            lam, mu = lam * Fraction(3, 4), mu * Fraction(3, 4) + P["L"] * ALPHA
        else:
            f = 1 - P["q"] * ALPHA
            lam, mu = lam * f, mu * f + P["c"] * P["q"] * ALPHA
    return lam, mu


def certify_word(K: int, x0: Fraction, w: str, xs=None, it=None,
                 i0: int = 0, delta_frac: Fraction = Fraction(1, 64)):
    """Interval-invariance certificate for a given word on the exact map.

    Returns dict(ok=..., phase, delta, dmin, lam, xstar, cycle, t_entry).
    `xs`/`it` (exact orbit) optional for entry-time reporting.
    """
    P = boundary_params(K)
    lam, mu = compose_word(w, P)
    if lam >= 1:
        return {"ok": False, "reason": "non-contracting word"}
    xstar = mu / (1 - lam)
    n = len(w)
    ps = [xstar]
    for br in w:
        ps.append(forced_step(ps[-1], P, br))
    if ps[n] != ps[0]:
        return {"ok": False, "reason": "fixed point does not realize word"}
    if any((p > P["h"]) != (br == "A") for p, br in zip(ps, w)):
        return {"ok": False, "reason": "cycle point on wrong side of h"}
    dmin = min(abs(p - P["h"]) for p in ps[:n])
    delta = delta_frac * dmin
    for j0 in range(n):
        lo, hi = ps[j0] - delta, ps[j0] + delta
        a, bb = lo, hi
        good = True
        for br in (w[j0:] + w[:j0]):
            if br == "A" and a > P["h"]:
                a, bb = (Fraction(3, 4) * a + P["L"] * ALPHA,
                         Fraction(3, 4) * bb + P["L"] * ALPHA)
            elif br == "B" and bb <= P["h"]:
                a, bb = (a + (P["c"] - a) * P["q"] * ALPHA,
                         bb + (P["c"] - bb) * P["q"] * ALPHA)
            else:
                good = False
                break
        if good and lo <= a and bb <= hi:
            t_entry = None
            if xs is not None:
                t = i0 + j0
                while t < len(xs):
                    if lo <= xs[t] <= hi:
                        t_entry = t
                        break
                    t += n
            return {"ok": True, "phase": j0, "delta": delta, "dmin": dmin,
                    "lam": lam, "xstar": xstar, "cycle": ps[:n],
                    "t_entry": t_entry}
    return {"ok": False, "reason": "no phase admitted interval",
            "dmin": dmin}


def full_certificate(K: int, x0: Fraction, T: int = 4000, rep: int = 20):
    """Detect (strongly) then certify. Returns dict with word data."""
    P = boundary_params(K)
    xs, it = scalar_orbit(x0, P, T)
    i0, n, w = detect_eventual_word(it, rep=rep)
    if w is None:
        return {"ok": False, "stage": "detection",
                "reason": f"no >= {rep}-repeat word within T={T}"}
    cert = certify_word(K, x0, w, xs=xs, it=it, i0=i0)
    cert.update({"t0": i0, "n": n, "word": w})
    return cert


# ---------------------------------------------------------------------------
# phase-map analysis (Lemma A/B) and pure windows (r = 3)
# ---------------------------------------------------------------------------

def pure_window_membership(K: int):
    """For r in {3,4}: which pure A^k B window does sigma(K) belong to?
    Window k: sigma in ((4/3)^(k-1) - (3/4) beta, (4/3)^k - beta].
    Returns k or None."""
    P = boundary_params(K)
    if K % 6 not in (3, 4):
        raise ValueError("pure windows defined for r in {3,4}")
    s, beta = P["sigma"], P["beta"]
    for k in range(1, 16):
        lo = Fraction(4, 3) ** (k - 1) - Fraction(3, 4) * beta
        hi = Fraction(4, 3) ** k - beta
        if lo < s <= hi:
            return k
    return None


def sigma_strictly_monotone(K_lo: int, K_hi: int, r: int) -> bool:
    """Check sigma(K) strictly monotone (direction inferred from first pair)
    on the class-r lattice between K_lo and K_hi (exact)."""
    Ks = [K for K in range(K_lo, K_hi + 1) if K % 6 == r]
    vals = [boundary_params(K)["sigma"] for K in Ks]
    dec = all(v2 < v1 for v1, v2 in zip(vals, vals[1:]))
    inc = all(v2 > v1 for v1, v2 in zip(vals, vals[1:]))
    return dec or inc
