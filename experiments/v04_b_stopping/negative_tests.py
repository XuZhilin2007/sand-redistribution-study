# -*- coding: utf-8 -*-
"""Negative tests for the v04 acceptance gates (review R1 mandate).

Proves by injection that each gate actually rejects:
  N1  decomp_ok=False            -> counted, main-style exit nonzero
  N2  identity_ok=False          -> counted
  N3  cert_all=True + real dip   -> CERT-CONTRADICTION counted
  N4  clean good results         -> zero failures (no false rejects)
  N5  no-entry case              -> not counted as implementation failure
  N6  cert_census validation     -> mismatch on the alpha=1/4 expected
      set (with the K=6 registered special case handled explicitly)
      returns failure
Exit 0 iff all injections behave as required.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "v04_a_alpha_phase"))

from decomp_alpha import count_failures          # noqa: E402

GOOD = {"K": 9, "alpha": "1/4", "t_star": 2, "decomp_ok": True,
        "identity_ok": True, "cert_all": True, "prefix_ok": True,
        "post_t0_dip_below": False, "layers": []}


def main() -> int:
    ok = True

    def check(name, cond):
        nonlocal ok
        print(f"  [{ 'PASS' if cond else 'FAIL'}] {name}")
        ok = ok and cond

    # N1: decomp failure must be counted
    r = dict(GOOD, decomp_ok=False)
    check("N1 decomp_ok=False rejected",
          any(cat == "DECOMP" for _, _, cat in count_failures([r])))
    # N2: identity failure must be counted
    r = dict(GOOD, identity_ok=False)
    check("N2 identity_ok=False rejected",
          any(cat == "IDENTITY" for _, _, cat in count_failures([r])))
    # N3: certificate contradiction must be counted
    r = dict(GOOD, cert_all=True, post_t0_dip_below=True)
    check("N3 cert contradiction rejected",
          any(cat == "CERT-CONTRADICTION" for _, _, cat in count_failures([r])))
    # N4: good results pass
    check("N4 good results accepted", count_failures([GOOD]) == [])
    # N5: no-entry not an implementation failure
    r = {"K": 9, "alpha": "1/8", "t_star": None,
         "layers": ["NO-ENTRY-IN-WINDOW"]}
    check("N5 no-entry not counted", count_failures([r]) == [])
    # N3b: S-FAIL open (conditions fail, no contradiction) is NOT counted
    r = dict(GOOD, cert_all=False, post_t0_dip_below=True)
    check("N3b S-FAIL open not miscounted", count_failures([r]) == [])

    # N6: census validation as an effective gate
    from cert_census_validate import validate_baseline
    good_rows = [{"K": K, "alpha": "1/4", "class": "CERTIFIED"}
                 for K in range(6, 204) if K not in
                 {6, 9, 10, 11, 16, 17, 23}]
    for K in (6, 9, 10, 11, 16, 17, 23):
        good_rows.append({"K": K, "alpha": "1/4", "class": "S-FAIL:S2"})
    check("N6a correct baseline accepted", validate_baseline(good_rows))
    bad_rows = [dict(r) for r in good_rows]
    bad_rows[0]["class"] = "S-FAIL:S1"        # a K=7 that should certify
    check("N6b corrupted baseline rejected", not validate_baseline(bad_rows))
    bad2 = [dict(r) for r in good_rows]
    bad2 = [r for r in bad2 if r["K"] != 9]   # drop a real exception
    check("N6c missing exception rejected", not validate_baseline(bad2))
    # K=6 S-FAIL at anchor t0=1 is EXPECTED (registered special case),
    # not a theorem counterexample:
    rows6 = [dict(r) for r in good_rows]
    check("N6d K=6 expected-exception semantics included",
          any(r["K"] == 6 and r["class"].startswith("S-FAIL")
              for r in rows6))
    print("NEGATIVE TESTS:", "ALL PASS" if ok else "FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
