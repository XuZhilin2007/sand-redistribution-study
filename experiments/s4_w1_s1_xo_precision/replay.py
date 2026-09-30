"""S1 high-precision replay engine (Stage 4 Wave 1 S1).

Re-implements the frozen M1A semantics (run_m1a_deterministic, adaptive.py)
in three arithmetic tiers so the XO / KC dynamics can be compared across
precision:

  * "exact" -- fractions.Fraction. Exact real-arithmetic trajectory. Legal
    because every bin probability involved is exactly rational at the
    execution grids: the canonical initial CDF is 2x - x^2 (rational at
    rational edges) and the XO respray CDF is piecewise linear with
    rational coefficients and rational breakpoints. The ppf radicals live
    only on the Monte-Carlo path, which the deterministic expected-mass
    model never touches. The stopping tolerance becomes the exact rational
    1/10^12; ties (worst == 0) are decidable exactly.
  * "mp50" / "mp100" -- mpmath at dps 50 / 100. Same mathematical
    semantics; bin probabilities are NOT renormalized (the model's
    telescoping sum is exactly 1; only the float64 implementation carries
    the fp residue-correction normalization).
  * float64 -- NOT re-implemented here: the F tier is the committed numpy
    implementation, regenerated through the committed runner modules in
    run_s1.py and string-verified against the committed CSVs (STAGE V).

Semantics (identical across tiers, matching run_m1a_deterministic):

    p_0 = initial bin probabilities
    each round t:  D_j = cumsum(p)_j - j/K   (j = 1..K)
                   s_t = min argmax D_j,  d_t = D_{s_t}
    stop iff d_t <= tol  (recorded at that t; no sweep executed)
    otherwise (t < H):  F_s = sum(p[:s]),  M = alpha * F_s
        p'_j = (1 - alpha) p_j + M q_j   (j <= s)
        p'_j = p_j + M q_j               (j >  s)

Diagnostics mirror the committed runner composition exactly
(kernel_theory.target_matched_excess / cumulative_mismatch /
rebound_residual; R_t over the index set {j : gamma_j > 0}, guard 0).

Exact identity used as a self-check (identity_check=True): by the M1C.2
one-step decomposition D'_j = D_TM,j + M*Delta_j, so

    max_j E_t(j) = max_j D'_{t+1,j} - d_t = d_{t+1} - d_t,

i.e. the worst_excess diagnostic IS the (exact) d-increment. The exact
tier verifies this identity per round with Fraction equality.
"""

from __future__ import annotations

import time
from fractions import Fraction
from typing import Sequence

import mpmath as mp

EXACT_TOL = Fraction(1, 10 ** 12)
STORAGE_DIGITS = 60
STORAGE_DPS = 80  # workspace for exact->decimal storage conversion
GUARD_DIGITS = 15  # mp tiers compute at requested dps + guard


# --------------------------------------------------------------------------
# scalar CDFs, per arithmetic
# --------------------------------------------------------------------------

def near_cdf_exact(x: Fraction) -> Fraction:
    """Canonical q_near CDF 2x - x^2 at an exact rational point."""
    return 2 * x - x * x


def xo_cdf_exact(x: Fraction) -> Fraction:
    """Frozen XO CDF at an exact rational point (piecewise linear,
    breakpoints 1/3, 2/3, 5/6; G_O = 1 on [5/6, 1])."""
    if x <= 0:
        return Fraction(0)
    if x < Fraction(1, 3):
        return x / 2
    if x <= Fraction(2, 3):
        return 3 * x / 2 - Fraction(1, 3)
    if x <= Fraction(5, 6):
        return 2 * x - Fraction(2, 3)
    return Fraction(1)


def near_cdf_mp(x):
    return 2 * x - x * x


def xo_cdf_mp(x):
    third = mp.mpf(1) / 3
    two_thirds = mp.mpf(2) / 3
    five_sixths = mp.mpf(5) / 6
    if x <= 0:
        return mp.mpf(0)
    if x < third:
        return x / 2
    if x <= two_thirds:
        return 3 * x / 2 - mp.mpf(1) / 3
    if x <= five_sixths:
        return 2 * x - two_thirds
    return mp.mpf(1)


CDF_TABLE = {
    "near": (near_cdf_exact, near_cdf_mp),
    "XO": (xo_cdf_exact, xo_cdf_mp),
}


# --------------------------------------------------------------------------
# bin probabilities
# --------------------------------------------------------------------------

def exact_bin_probs(cdf, K: int) -> list[Fraction]:
    """q_j = G(j/K) - G((j-1)/K) in exact rational arithmetic.

    Asserts the identities that hold exactly for a CDF: q_j >= 0 and the
    telescoping sum == 1 (this replaces the float implementation's residue
    normalization, which has no exact content).
    """
    edges = [Fraction(j, K) for j in range(K + 1)]
    q = [cdf(edges[j]) - cdf(edges[j - 1]) for j in range(1, K + 1)]
    if any(qj < 0 for qj in q):
        raise ValueError("exact bin probs: negative bin mass")
    if sum(q) != 1:
        raise ValueError("exact bin probs: telescoping sum != 1")
    return q


def mp_bin_probs(cdf, K: int, dps: int) -> list:
    """mp tier bin probabilities: CDF differences on exact rational edges
    converted at full workspace precision. No renormalization."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        q = []
        prev = mp.mpf(0)
        for j in range(1, K + 1):
            cur = cdf(mp.mpf(j) / mp.mpf(K))
            q.append(cur - prev)
            prev = cur
        if any(qj < 0 for qj in q):
            raise ValueError("mp bin probs: negative bin mass")
        return q
    finally:
        mp.mp.dps = old


def scalar_bin_probs(kernel_law: str, K: int, tier: str, work_dps: int = 0) -> list:
    if tier == "exact":
        return exact_bin_probs(CDF_TABLE[kernel_law][0], K)
    if tier not in ("mp50", "mp100"):
        raise ValueError(f"unknown tier {tier!r}")
    return mp_bin_probs(CDF_TABLE[kernel_law][1], K, work_dps)


# --------------------------------------------------------------------------
# one-round semantics + diagnostics (tier-polymorphic; one = tier's 1)
# --------------------------------------------------------------------------

def cumsum_excess(p: Sequence, K: int) -> list:
    """D_j = cumsum(p)_j - j/K in the tier's own arithmetic."""
    zero = p[0] - p[0]
    one = zero + 1
    D = []
    acc = zero
    for j, pj in enumerate(p, start=1):
        acc = acc + pj
        D.append(acc - one * j / K)
    return D


def argmax_first(D: Sequence) -> int:
    """1-based index of the first maximal entry (frozen min-argmax rule)."""
    best = 0
    for j in range(1, len(D)):
        if D[j] > D[best]:
            best = j
    return best + 1


def update_state(p: Sequence, q: Sequence, s: int, alpha) -> list:
    """One active sweep: prefix scaled by (1-alpha), respray M*q everywhere,
    M = alpha * F_s from the pre-update state."""
    zero = p[0] - p[0]
    F_s = zero
    for j in range(s):
        F_s = F_s + p[j]
    M = alpha * F_s
    one_minus = (zero + 1) - alpha
    return [(one_minus * pj if j < s else pj) + M * qj
            for j, (pj, qj) in enumerate(zip(p, q))]


def diagnostics(p: Sequence, q: Sequence, K: int, s: int, alpha) -> dict:
    """worst_t = max_j E_t(j) and R_t, mirroring the committed composition:
    E = M * Delta - gamma with gamma = d_max - D_TM-profile;
    R_t = max over {j: gamma_j > 0} of M * max(Delta_j, 0) / gamma_j."""
    zero = p[0] - p[0]
    one = zero + 1
    F = []
    acc = zero
    for pj in p:
        acc = acc + pj
        F.append(acc)
    F_s = F[s - 1]
    M = alpha * F_s
    T = [one * j / K for j in range(1, K + 1)]
    delta = []
    acc = zero
    for qj in q:
        acc = acc + qj
        delta.append(acc)
    delta = [delta[j] - T[j] for j in range(K)]
    D = [F[j] - T[j] for j in range(K)]
    d_max = D[s - 1]
    D_tm = [((one - alpha) * D[j] - alpha * T[j] * (one - F_s)) if j < s
            else (D[j] - M * (one - T[j])) for j in range(K)]
    gamma = [d_max - v for v in D_tm]
    E = [M * delta[j] - gamma[j] for j in range(K)]
    worst = E[0]
    for j in range(1, K):
        if E[j] > worst:
            worst = E[j]
    R = zero
    for j in range(K):
        if gamma[j] > 0:
            dpos = delta[j] if delta[j] > 0 else zero
            r = M * dpos / gamma[j]
            if r > R:
                R = r
    return {"worst": worst, "R": R, "moved": M}


# --------------------------------------------------------------------------
# full replay
# --------------------------------------------------------------------------

def anchor_m_b(K: int) -> tuple[int, Fraction]:
    """(m, b = m/K) with m = ceil(5K/6); asserts the ceiling identity."""
    m = -(-5 * K // 6)
    if not (5 * K <= 6 * m < 5 * K + 6):
        raise AssertionError("ceil(5K/6) formula broken")
    return m, Fraction(m, K)


def replay(kernel: str, K: int, tier: str, H: int,
           identity_check: bool = False) -> dict:
    """Full trajectory replay in the requested tier.

    kernel: "XO" (canonical init, XO respray) or "KC" (canonical both).
    tier:   "exact" | "mp50" | "mp100".  Horizon H = number of active
    rounds; a stop records stopped_at and stops before sweeping.
    Rows carry raw tier values under underscore keys (never written to
    CSV) plus storage-precision strings for output.
    """
    t0 = time.perf_counter()
    if kernel not in ("XO", "KC"):
        raise ValueError(kernel)
    if tier == "exact":
        alpha, tol, boundary_tol = Fraction(1, 4), EXACT_TOL, EXACT_TOL
        mp.mp.dps = STORAGE_DPS  # storage-conversion workspace only
    elif tier in ("mp50", "mp100"):
        # requested precision + guard digits; storage stays within the
        # requested part (60 digits < 50/100)
        mp.mp.dps = (50 if tier == "mp50" else 100) + GUARD_DIGITS
        alpha, tol, boundary_tol = mp.mpf(1) / 4, mp.mpf(10) ** -12, mp.mpf(10) ** -12
    else:
        raise ValueError(tier)
    init_law, throw_law = ("near", "XO") if kernel == "XO" else ("near", "near")
    p = scalar_bin_probs(init_law, K, tier, mp.mp.dps)
    q = scalar_bin_probs(throw_law, K, tier, mp.mp.dps)
    m, b = anchor_m_b(K)
    anchor = b * (1 - b)
    if tier != "exact":
        # never hand a Fraction to mpmath directly: its conversion goes
        # through float64; convert at the working precision instead
        anchor = mp.mpf(anchor.numerator) / mp.mpf(anchor.denominator)
    rows: list[dict] = []
    stopped_at = None
    exact_flags_ok = True
    identity_ok = True

    def check_identity(prev_row: dict, D_next: Sequence) -> None:
        nonlocal identity_ok
        if prev_row is None:
            return
        d_next_max = D_next[argmax_first(D_next) - 1]
        if tier == "exact":
            ok = d_next_max - prev_row["_d"] == prev_row["_worst"]
        else:
            ok = abs(d_next_max - prev_row["_d"] - prev_row["_worst"]) \
                <= mp.mpf(10) ** (-(mp.mp.dps - 10))
        identity_ok = identity_ok and ok

    prev_row = None
    for t in range(H):
        D = cumsum_excess(p, K)
        if identity_check:
            check_identity(prev_row, D)
        s = argmax_first(D)
        d = D[s - 1]
        if d <= tol:
            stopped_at = t
            break
        if tier == "exact" and not d > EXACT_TOL:
            raise AssertionError("exact tier: stop-rule sanity broken")
        dg = diagnostics(p, q, K, s, alpha)
        sbound_ok = s <= m
        inv_ok = True
        if tier == "exact":
            inv_ok = D[m - 1] == anchor
            exact_flags_ok = exact_flags_ok and inv_ok and sbound_ok
        row = {
            "t": t, "s": s,
            "d": _to_storage(d), "e": _to_storage(d - anchor),
            "moved": _to_storage(dg["moved"]),
            "worst": _to_storage(dg["worst"]), "R": _to_storage(dg["R"]),
            "boundary": abs(dg["worst"]) <= boundary_tol,
            "tie_exact": (dg["worst"] == 0) if tier == "exact" else None,
            "invariant_exact_ok": inv_ok if tier == "exact" else None,
            "s_bound_ok": sbound_ok,
            "_d": d, "_worst": dg["worst"],
        }
        rows.append(row)
        prev_row = row
        p = update_state(p, q, s, alpha)
    if identity_check and rows:
        D_final = cumsum_excess(p, K)
        check_identity(prev_row, D_final)
    return {
        "rows": rows, "stopped_at": stopped_at, "censored": stopped_at is None,
        "K": K, "kernel": kernel, "tier": tier, "m": m, "anchor": anchor,
        "exact_flags_ok": exact_flags_ok, "identity_ok": identity_ok,
        "wall_seconds": round(time.perf_counter() - t0, 2),
        "final_state": p,
    }


def rebound_flags(result: dict) -> list[bool]:
    """Strict d_{t+1} > d_t flags over all active rounds (no tolerance),
    matching the committed runners' post-hoc strict comparison."""
    rows = result["rows"]
    flags = [b["_d"] > a["_d"] for a, b in zip(rows, rows[1:])]
    if rows:
        D_final = cumsum_excess(result["final_state"], result["K"])
        flags.append(D_final[argmax_first(D_final) - 1] > rows[-1]["_d"])
    return flags


def _to_storage(x):
    """~60 significant digits for CSV storage. Fraction -> mp conversion
    must go through numerator/denominator at a raised workspace (mpmath's
    direct Fraction conversion passes through float64)."""
    if isinstance(x, Fraction):
        old = mp.mp.dps
        mp.mp.dps = STORAGE_DPS
        try:
            x = mp.mpf(x.numerator) / mp.mpf(x.denominator)
        finally:
            mp.mp.dps = old
    return mp.nstr(x, STORAGE_DIGITS)
