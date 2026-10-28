# -*- coding: utf-8 -*-
"""Coverage/STRUCTURE validation for interval-certificate result sets.

SCOPE (review round-3 clarification): this module does NOT re-verify
any mathematical inequality.  It only checks that records cover the
caller-supplied target K set exactly once, carry the expected engine
version, valid classes and required fields, and that interval bounds,
modes and time windows are structurally legal (round-5 freeze guard:
0 < lo < hi < 1 with parseable rational endpoints, defined modes,
T1 a positive integer; tie-split records validated per their real
format with legal partial fields allowed for OPEN).  Mathematical
acceptance is a separate layer -- the independent auditor
tools/audit_v04_interval_v3.py (exact Bernstein coefficient checks of
every sign obligation).  A record set passing this validator is NOT
thereby mathematically certified.

Returns list of problems (empty == structurally valid)."""
from __future__ import annotations

from fractions import Fraction

VALID_CLASSES = {"PASS", "OPEN"}
VALID_MODES = {"closed", "leftopen", "tie-split"}
REQUIRED_FIELDS = ("engine", "K", "lo", "hi", "mode", "T1", "class",
                   "fails", "coverage")
# single source of truth for the PRODUCTION engine version (round-3
# fix 3.1: interval_certs imports this; legacy v3 records are accepted
# ONLY through the independent mathematical auditor with an explicit
# engine="v04-interval-v3" structural pass, never by editing fields)
PRODUCTION_ENGINE = "v04-interval-v4"
CURRENT_ENGINE = PRODUCTION_ENGINE


def _parse_frac(s):
    try:
        return Fraction(s)
    except (ValueError, TypeError, ZeroDivisionError):
        return None


def _bad_t1(t1):
    # closed/leftopen: positive int; tie-split: [pos int, pos int] or
    # [pos int, None]/[None, pos int] for OPEN partials are allowed --
    # but PASS requires both present.
    return False


def _check_common(rec, problems, engine):
    K = rec.get("K")
    for f in REQUIRED_FIELDS:
        if f not in rec:
            problems.append(f"K={K}: missing field {f}")
    if rec.get("engine") != engine:
        problems.append(f"K={K}: engine {rec.get('engine')} != {engine}")
    if rec.get("class") not in VALID_CLASSES:
        problems.append(f"K={K}: bad class {rec.get('class')}")


def _check_piece(rec, problems):
    K = rec.get("K")
    lo = _parse_frac(rec.get("lo"))
    hi = _parse_frac(rec.get("hi"))
    if lo is None:
        problems.append(f"K={K}: unparseable lo {rec.get('lo')!r}")
        lo = None
    if hi is None:
        problems.append(f"K={K}: unparseable hi {rec.get('hi')!r}")
    if lo is not None and hi is not None:
        if not (0 < lo < hi < 1):
            problems.append(f"K={K}: reversed/illegal bounds "
                            f"[{rec.get('lo')},{rec.get('hi')}]")
    mode = rec.get("mode")
    if mode not in VALID_MODES - {"tie-split"}:
        problems.append(f"K={K}: undefined mode {mode!r}")
    t1 = rec.get("T1")
    if not (isinstance(t1, int) and t1 > 0):
        problems.append(f"K={K}: illegal T1 {t1!r}")


def _check_tiesplit(rec, problems):
    K = rec.get("K")
    if rec.get("mode") != "tie-split":
        return
    # real tie-split format: lo/hi legal fractions, T1 a 2-list of
    # positive ints (None allowed only in OPEN records), fails a dict
    # with left/right keys
    lo = _parse_frac(rec.get("lo"))
    hi = _parse_frac(rec.get("hi"))
    if lo is None or hi is None or not (0 < lo < hi < 1):
        problems.append(f"K={K}: tie-split illegal bounds")
    t1 = rec.get("T1")
    if not (isinstance(t1, list) and len(t1) == 2):
        problems.append(f"K={K}: tie-split T1 not a 2-list: {t1!r}")
    else:
        for x in t1:
            if x is None and rec.get("class") == "OPEN":
                continue
            if not (isinstance(x, int) and x > 0):
                problems.append(f"K={K}: tie-split illegal T1 entry {x!r}")
    if not isinstance(rec.get("fails"), dict):
        problems.append(f"K={K}: tie-split fails not a dict")


def validate_interval_rows(rows, target_KS, engine=CURRENT_ENGINE):
    problems = []
    counts = {}
    for rec in rows:
        K = rec.get("K")
        counts[K] = counts.get(K, 0) + 1
        _check_common(rec, problems, engine)
        if rec.get("mode") == "tie-split":
            _check_tiesplit(rec, problems)
        else:
            _check_piece(rec, problems)
    for K in target_KS:
        if counts.get(K, 0) == 0:
            problems.append(f"target K={K} missing")
        elif counts.get(K, 0) > 1:
            problems.append(f"target K={K} duplicated")
    for K in counts:
        if K not in set(target_KS):
            problems.append(f"non-target K={K} present")
    return problems
