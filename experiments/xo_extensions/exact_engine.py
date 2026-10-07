"""Portable exact checks for the canonical XO supplement.

The mass-vector engine is carried over from the recorded independent review.
All H choices are pre-update states t=0,...,H-1. No scalar reduction advances
the integer engine. These finite checks do not prove universal quantifiers.
"""
from __future__ import annotations

import ast
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FROZEN_S2 = "experiments/s4_w1_s2_boundary_microscope/run_s2.py"


def frozen_classifier(root=None):
    """Load the deposited classifier without importing/running the S2 CLI."""
    root = ROOT if root is None else Path(root)
    source = root / FROZEN_S2
    tree = ast.parse(source.read_text(encoding="utf-8"))
    node = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "classify_itinerary")
    scope = {}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[])),
                 FROZEN_S2, "exec"), scope)
    return scope["classify_itinerary"]


def boundary_step(K, alpha, x):
    """Exact boundary observable update; valid after boundary confinement."""
    assert K >= 6 and 0 < alpha < 1
    m = (5 * K + 5) // 6
    h, b = F(1, K), F(m, K)
    c, a, q = 2 * b - b * b, b * (1 - b), F(6 - K % 6, 3 * K)
    if x > h:
        return m, a, (1 - alpha) * x + alpha * c * q
    return m - 1, a + h - x, x + alpha * (c - x) * q



def cdf_num(j, K):
    """G(j/K) with common denominator 6K, using integer comparisons."""
    if 3 * j <= K:
        return 3 * j
    if 3 * j <= 2 * K:
        return 9 * j - 2 * K
    if 6 * j <= 5 * K:
        return 12 * j - 4 * K
    return 6 * K


def integer_orbit(K, alpha, H, frozen):
    """Evolve all mass numerators with a shared integer denominator.

    This uses neither the lab's Fraction-vector engine nor its CDF engine.
    No boundary reduction is used to advance the full trajectory.
    """
    m = (5 * K + 5) // 6
    den = K * K
    nums = [2 * K - 2 * j + 1 for j in range(1, K + 1)]
    gden = 6 * K
    gn = [cdf_num(j, K) - cdf_num(j - 1, K) for j in range(1, K + 1)]
    anum, aden = alpha.numerator, alpha.denominator
    h, b = F(1, K), F(m, K)
    c, a, q = 2 * b - b * b, b * (1 - b), F(gn[m - 1], gden)
    itinerary, first, confinement = [], None, None
    boundary_equalities = 0
    closure_states = 0
    for t in range(H):
        assert min(nums) >= 0 and sum(nums) == den
        pref, prefs, scores = 0, [], []
        for j, num in enumerate(nums, 1):
            pref += num
            prefs.append(pref)
            scores.append(K * pref - j * den)
        largest = max(scores)
        s = scores.index(largest) + 1
        itinerary.append(s)
        assert s <= m and prefs[m - 1] * c.denominator == c.numerator * den
        assert largest * a.denominator >= a.numerator * K * den
        assert all(z * a.denominator < a.numerator * K * den for z in scores[m:])
        internal_below = all(z * a.denominator < a.numerator * K * den for z in scores[:m - 2])
        if internal_below and confinement is None:
            confinement = t
        if s in (m - 1, m) and first is None:
            first = t
        if first is not None:
            assert s in (m - 1, m), (K, str(alpha), t, "interior return")
            assert s == (m if K * nums[m - 1] > den else m - 1)
            discrepancy_excess = max(den - K * nums[m - 1], 0)
            assert (largest - discrepancy_excess) * a.denominator == a.numerator * K * den
            boundary_equalities += K * nums[m - 1] == den
            closure_states += 1
        moved_num = anum * prefs[s - 1]
        new_nums = [((aden - anum) if j < s else aden) * n * gden + moved_num * gn[j]
                    for j, n in enumerate(nums)]
        new_den = den * aden * gden
        if first is not None:
            x = F(nums[m - 1], den)
            predicted = (1 - alpha) * x + alpha * c * q if x > h else x + alpha * (c - x) * q
            assert new_nums[m - 1] * predicted.denominator == predicted.numerator * new_den
        nums, den = new_nums, new_den
    bound = int(F(81 * K, 80) / alpha) + 1
    if first is not None:
        assert first <= bound
    if confinement is not None:
        assert confinement <= bound
    value = frozen(itinerary, m, K)
    return {"K": K, "alpha": str(alpha), "H": H, **value, "confinement_time": confinement,
            "entry_bound": bound, "boundary_equalities": boundary_equalities,
            "closure_states": closure_states,
            "itinerary_sha256": hashlib.sha256(json.dumps(itinerary).encode()).hexdigest()}, itinerary


def proof_identities():
    records = []
    for K in list(range(6, 401)) + [1000, 1001, 1002, 1003, 1004, 1005, 10**12 + 3]:
        h, r = F(1, K), K % 6
        m = (5 * K + 5) // 6
        b = F(m, K)
        c, a, q = 2 * b - b * b, b * (1 - b), F(6 - r, 3 * K)
        assert F(5, 6) <= b <= F(10, 11) and a >= F(10, 121) and c >= F(35, 36)
        assert q == F(cdf_num(m, K) - cdf_num(m - 1, K), 6 * K)
        # Check all interior coordinates for modest K, and segment endpoints for huge K.
        points = range(1, m - 1) if K <= 1005 else sorted({1, K // 3, 2 * K // 3, m - 2})
        for j in points:
            y = j * h
            gj = F(cdf_num(j, K), 6 * K)
            assert a - (c * gj - y) >= F(29, 108) * h
            assert 1 - gj >= 2 * h
        assert min(F(29, 108) * h, 2 * a * h) >= F(20, 121) / K
        if K >= 7:
            z = (c - h) * (1 - q - h / 2) - b + 2 * h
            u = K // 6 - 1
            polynomial = 180 * u**3 + (231 + 136 * r) * u**2 + (12 + 152 * r + 30 * r**2) * u + 2 * r**3 + 19 * r**2 + 31 * r - 39
            assert 6 * K**3 * z == polynomial and polynomial > 0
        if r < 3:
            assert c * q > h and F(2011, 580) * K - 2 >= K
        else:
            assert c * q / h <= F(80, 81)
            assert 81 * K * F(3, 4)**(3 * K) < 1 if K <= 1005 else True
            assert F(72, 29) * K + 1 <= 3 * K
            assert F(141, 20) * K + 12 < 10 * K - 1
            assert (h - c * q) / (h - F(3, 4) * c * q) >= F(1, 21)
        records.append(K)
    assert F(3, 4)**11 <= F(1, 21)
    assert F(10, 9) * F(27, 64) == F(15, 32)
    return {"status": "PASS", "K_count": len(records), "K_values": records,
            "coverage_note": "Finite algebra checks support the separately reviewed general proofs; no universal quantifier is inferred from this scan."}
