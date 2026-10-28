# -*- coding: utf-8 -*-
"""Negative/boundary regressions for the v3 interval certifier
(review round-2 mandates): I1 (symkernel selftest), I2 end-to-end
wrong-proposition rejection, W0 necessity, coverage validator I6-style,
Gram identity on real assemblies.
Exit 0 iff all behave as required."""
from __future__ import annotations

import sys
import json
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import interval_certs as ic          # noqa: E402
import symkernel                     # noqa: E402
from validate_interval import validate_interval_rows  # noqa: E402


def true_U1_minus_U2(K: int, alpha: F):
    """Independent exact replay value (small K, direct)."""
    def G(x):
        if x < F(1, 3):
            return x / 2
        if x <= F(2, 3):
            return 3 * x / 2 - F(1, 3)
        if x <= F(5, 6):
            return 2 * x - F(2, 3)
        return F(1)
    p = [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]
    g = [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]
    states = [p]
    for _ in range(2):
        cur = states[-1]
        acc, best, Dbest = F(0), 0, None
        for j, pj in enumerate(cur):
            acc += pj
            d = acc - F(j + 1, K)
            if Dbest is None or d > Dbest:
                Dbest, best = d, j
        s = best + 1
        M = alpha * sum(cur[:s])
        states.append([((1 - alpha) * pj if j < s else pj) + M * gj
                       for j, (pj, gj) in enumerate(zip(cur, g))])
    U = [K * sum((y - F(1, K)) ** 2 for y in st) for st in states]
    return U[1] - U[2]


def main() -> int:
    ok = True

    def check(name, cond):
        nonlocal ok
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
        ok = ok and cond

    # I1: symbol kernel boundary tests (21 checks, incl. review poly)
    check("I1 symkernel selftest", symkernel.sturm_selftest())

    # I2: K=6, r=1/20 -> interval [1/5, 3/10], T1=5 must be OPEN with
    # the W0 prefix failure (was wrongly PASS in v2)
    val = true_U1_minus_U2(6, F(3, 10))
    check("I2: true U1-U2 at a=3/10 equals -11669/2160000",
          val == F(-11669, 2160000))
    rec = ic.analyze_piece(6, F(1, 5), F(3, 10), "closed", F(1, 4), 5)
    check("I2: analyze_piece(6,[1/5,3/10]) is OPEN",
          rec["class"] == "OPEN")
    check("I2: failure is the W0 prefix (not a spurious cause)",
          "W0:U1<=U2" in rec["fails"])

    # W0 also guards smaller violations: K=6 at r=1/8 must be OPEN
    rec = ic.analyze_piece(6, F(1, 8), F(3, 8), "closed", F(1, 4), 8)
    check("W0 guards: K=6 r=1/8 OPEN", rec["class"] == "OPEN")

    # a genuinely valid small case still passes (no over-rejection)
    rec = ic.certify_symmetric(6, F(1, 125), [15])
    check("K=6 r=1/125 still certifies", rec["class"] == "PASS")
    rec = ic.certify_symmetric(9, F(1, 500), [15])
    check("K=9 r=1/500 still certifies (E-exception member)",
          rec["class"] == "PASS")

    # I5: Gram identity present on real assemblies
    check("I5 route tag", rec["route"] == "H-convex(Gram)")

    # I6-style coverage validator
    good = [{"K": K, "engine": ic.VERSION, "class": "PASS",
             "lo": "3/20", "hi": "1/5", "mode": "closed", "T1": 15,
             "fails": [], "word_head": [], "route": "H",
             "coverage": "x"} for K in ic.KS]
    check("coverage: correct rows accepted",
          validate_interval_rows(good, ic.KS) == [])
    miss = [r for r in good if r["K"] != 7]
    check("coverage: missing ordinary pass K rejected",
          len(validate_interval_rows(miss, ic.KS)) > 0)
    miss2 = [r for r in good if r["K"] != 1001]
    check("coverage: missing anchor rejected",
          len(validate_interval_rows(miss2, ic.KS)) > 0)
    dup = good + [dict(good[0])]
    check("coverage: duplicate rejected",
          len(validate_interval_rows(dup, ic.KS)) > 0)
    extra = good + [{"K": 204, "engine": ic.VERSION, "class": "PASS",
                     "lo": "3/20", "hi": "1/5", "mode": "closed",
                     "T1": 15, "fails": [], "word_head": [],
                     "route": "H", "coverage": "x"}]
    check("coverage: extra K rejected",
          len(validate_interval_rows(extra, ic.KS)) > 0)
    badv = [dict(r) for r in good]
    badv[3]["engine"] = "v2"
    check("coverage: stale engine version rejected",
          len(validate_interval_rows(badv, ic.KS)) > 0)

    # tie-split smoke: well-formed record for K=16
    rec = ic.certify_tie_split(16, F(1, 500), [15])
    check("tie-split K=16 well-formed",
          rec["mode"] == "tie-split" and "left" in rec["fails"])

    # ---- round-3 fixes: budget isolation, propagation, resume gate ----
    # 3.2a: consecutive top-level weak-nonneg calls do NOT share budget
    # (reviewer's case: tangent (y-1/4)^2 on [0, 1/2], 2100 calls)
    tangent = [F(1, 16), F(-1, 2), F(1)]
    okk = True
    for _ in range(2100):
        if not symkernel.nonneg_closed(tangent, F(0), F(1, 2)):
            okk = False
            break
    check("R3.2a no cross-call budget depletion (2100 calls)", okk)
    # 3.2b: BudgetExceeded in a piece returns OPEN (not a crash, not PASS)
    rec = ic.certify_symmetric(9, F(1, 500), [15])
    check("R3.2b certify returns well-formed record", "class" in rec)
    # 3.1: stale-engine record rejected by default; explicit legacy pass
    row_v3 = dict(good[0], engine="v04-interval-v3")
    check("R3.1 stale engine rejected by default validator",
          len(validate_interval_rows([row_v3], [row_v3["K"]])) > 0)
    check("R3.1 explicit legacy engine pass allowed",
          validate_interval_rows([row_v3], [row_v3["K"]],
                                 engine="v04-interval-v3") == [])
    check("R3.1 production validator accepts production records",
          validate_interval_rows(good, ic.KS) == [])
    # 3.3: bad resume record (missing field, BROKEN engine, reversed
    # bounds) must be structurally rejected
    bad_rec = {"K": 6, "engine": "BROKEN", "class": "PASS",
               "hi": "1/5", "lo": "3/10", "mode": "closed", "T1": 15,
               "fails": [], "word_head": [], "route": "x"}
    check("R3.3 malformed resume record rejected",
          len(validate_interval_rows([bad_rec], [6])) > 0)

    # ---- round-5 freeze guards: one-field-at-a-time mutations ------
    base = {"K": 6, "engine": ic.VERSION, "class": "PASS", "lo": "3/20",
            "hi": "1/5", "mode": "closed", "T1": 15, "fails": [],
            "word_head": [], "route": "H", "coverage": "x"}
    check("R5 legal record accepted",
          validate_interval_rows([base], [6]) == [])
    m1 = dict(base, lo="3/10")            # lo > hi
    check("R5 reversed bounds rejected",
          len(validate_interval_rows([m1], [6])) > 0)
    m2 = dict(base, lo="garbage")         # unparseable endpoint
    check("R5 unparseable lo rejected",
          len(validate_interval_rows([m2], [6])) > 0)
    m3 = dict(base, mode="garbage")       # undefined mode
    check("R5 undefined mode rejected",
          len(validate_interval_rows([m3], [6])) > 0)
    m4 = dict(base, T1=-10)               # illegal time window
    check("R5 negative T1 rejected",
          len(validate_interval_rows([m4], [6])) > 0)
    m5 = dict(base, T1=0)
    check("R5 zero T1 rejected",
          len(validate_interval_rows([m5], [6])) > 0)
    # tie-split real format: legal PASS + partial OPEN tolerated
    ts = {"K": 16, "engine": ic.VERSION, "class": "PASS", "lo": "61/250",
          "hi": "64/250", "mode": "tie-split", "T1": [15, 15],
          "fails": {"left": [], "right": []}, "word_head": {},
          "route": "lr", "coverage": "x"}
    check("R5 legal tie-split accepted",
          validate_interval_rows([ts], [16]) == [])
    ts_bad = dict(ts, T1=[15, -3])
    check("R5 tie-split bad T1 entry rejected",
          len(validate_interval_rows([ts_bad], [16])) > 0)

    # ---- round-5: REAL BudgetExceeded injection -> OPEN ------------
    import symkernel as sk

    def raise_budget(*a, **k):
        raise sk.BudgetExceeded()
    orig_gap = sk._gap_probes
    sk._gap_probes = raise_budget
    try:
        rec = ic.certify_symmetric(9, F(1, 500), [15])
        ok_inj = (rec["class"] == "OPEN"
                  and "BUDGET-EXCEEDED" in rec["fails"])
        rec_t = ic.certify_tie_split(16, F(1, 500), [15])
        ok_inj_t = (rec_t["class"] == "OPEN"
                    and any("BUDGET-EXCEEDED" in (v or [])
                            for v in (rec_t["fails"].get("left"),
                                      rec_t["fails"].get("right"))))
    finally:
        sk._gap_probes = orig_gap
    check("R5 budget injection -> symmetric OPEN w/ reason", ok_inj)
    check("R5 budget injection -> tie-split OPEN keeps reasons", ok_inj_t)

    # ---- round-5: bad records actually excluded from resume set ----
    from interval_certs import load_checkpoint
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        ck = Path(td) / "ck.jsonl"
        lines = [json.dumps(base), json.dumps(dict(base, lo="3/10")),
                 json.dumps(dict(base, T1=-10)),
                 json.dumps(dict(base, engine="v04-interval-v3"))]
        ck.write_text("\n".join(lines), encoding="utf-8")
        seen, rejected = load_checkpoint(ck, {6})
        check("R5 resume adopts ONLY the legal record",
              list(seen) == [6] and len(rejected) == 3)

    print("NEGATIVE TESTS (v3):", "ALL PASS" if ok else "FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
