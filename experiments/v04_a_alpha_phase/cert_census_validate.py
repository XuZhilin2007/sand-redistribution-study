# -*- coding: utf-8 -*-
"""Baseline validation for the v04 certificate census (review R1/R7).

At alpha=1/4 the census must reproduce the v0.3.0 certification
structure: S-FAIL exactly on {9,10,11,16,17,23} (theorem E finite
exception set) plus K=6 (registered special case: entry t*=1, covered
by its own point certificate -- an EXPECTED S-FAIL at anchor t0=1, not
a theorem counterexample); every other K CERTIFIED with no residence
flag.  Any mismatch is a validation FAILURE (used as an exit gate).
"""
from __future__ import annotations

from fractions import Fraction as F

BASE_ALPHA = F(1, 4)
EXPECTED_S_FAIL = {9, 10, 11, 16, 17, 23}   # E finite exception set
EXPECTED_SPECIAL = {6}                      # t*=1 registered special case


def validate_baseline(rows) -> bool:
    """rows: iterable of dicts with keys K, alpha (str), class."""
    sel = [r for r in rows if F(r["alpha"]) == BASE_ALPHA]
    if not sel:
        return False
    s_fail = {r["K"] for r in sel if r["class"].startswith("S-FAIL")}
    certified = {r["K"] for r in sel
                 if r["class"].startswith("CERTIFIED")
                 and "RES-FLAG" not in r["class"]}
    ks = {r["K"] for r in sel}
    return (s_fail == EXPECTED_S_FAIL | EXPECTED_SPECIAL
            and certified == ks - s_fail)
