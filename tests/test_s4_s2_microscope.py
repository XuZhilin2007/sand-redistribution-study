"""S4 Wave 1 / S2 microscope tests (frozen config e706a0d).

Run:  python tests/test_s4_s2_microscope.py

T1  reduction unit checks at K in {6, 7, 8, 9, 12} (short exact slices):
    branch equations (C3, measured F), selection/branch consistency (C4),
    d-prediction (C5) hold as exact Fraction identities post-entry;
    F_m conservation (C1) and the q_K closed form (C7) hold exactly;
T2  q_K closed form (6 - (K mod 6))/(3K) == definitional G_O(m/K)
    - G_O((m-1)/K) for ALL frozen K (census 6..53 + anchors), exact;
T3  frozen classifier unit tests on synthetic itineraries covering LOCK,
    CYCLE, and every OTHER subtype (incl. short_tail_lock, m_minus_1_only,
    interior_return, never_enters);
T4  Step-0 committed-CSV gate reproduces the freeze-time dry run
    (anchor classes LOCK/CYCLE/LOCK/CYCLE, all S0 checks pass).
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "s4_w1_s1_xo_precision"))
sys.path.insert(0, str(REPO / "experiments" / "s4_w1_s2_boundary_microscope"))

import numpy  # noqa: E402,F401  (engine import path sanity)

import replay  # noqa: E402
from run_s2 import (  # noqa: E402
    ALL_K,
    classify_itinerary,
    exact_replay_s2,
    step0,
)

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    if ok:
        print(f"[PASS] {name}" + (f" — {detail}" if detail else ""))
    else:
        print(f"[FAIL] {name} — {detail}")
        FAILURES.append(name)


# ---------------------------------------------------------------- T1
for K in (6, 7, 8, 9, 12):
    res = exact_replay_s2(K, instrument_full=False, H=25)
    m = res["m"]
    post = [r for r in res["raw"][res["cls"]["t_star"]:] ]
    c3 = all(r["branch_eq_ok"] for r in post)
    c4 = all(r["sel_branch_ok"] for r in post)
    c5 = all(r["dpred_ok"] for r in post)
    c1 = not res["c1_bad_states"]
    q_def = replay.exact_bin_probs(replay.xo_cdf_exact, K)[m - 1]
    c7 = (res["qK"] == q_def == Fraction(6 - (K % 6), 3 * K))
    check(f"T1 K={K}: C3 branch equations exact post-entry (slice H=25)", c3)
    check(f"T1 K={K}: C4 selection/branch consistency", c4)
    check(f"T1 K={K}: C5 d-prediction exact", c5)
    check(f"T1 K={K}: C1 F_m conservation over slice", c1)
    check(f"T1 K={K}: C7 q_K definitional == closed form", c7)

# ---------------------------------------------------------------- T2
ok_all = True
bad = []
for K in ALL_K:
    m, _b = replay.anchor_m_b(K)
    q_def = replay.exact_bin_probs(replay.xo_cdf_exact, K)[m - 1]
    q_closed = Fraction(6 - (K % 6), 3 * K)
    if not (q_def == q_closed):
        ok_all = False
        bad.append(K)
check("T2 q_K closed form == definitional for all 51 frozen K", ok_all,
      f"mismatches: {bad}")

# ---------------------------------------------------------------- T3
K, m = 6, 5
cases = [
    ([3, 4, 5, 5, 5, 5, 5, 5], "LOCK", ""),
    ([3, 4, 5, 4, 5, 4, 5, 4], "CYCLE", ""),
    ([3, 3, 3, 3, 3, 3], "OTHER", "never_enters"),
    ([3, 4, 3, 4, 5, 5, 5, 5], "OTHER", "interior_return"),
    ([3, 4, 4, 4, 4, 4, 4, 4], "OTHER", "m_minus_1_only"),
    ([3, 5, 5, 5], "OTHER", "short_tail_lock"),
]
for it, want_cls, want_sub in cases:
    got = classify_itinerary(it, m, K)
    check(f"T3 classify {it[:4]}... -> {want_cls}/{want_sub or '-'}",
          got["class"] == want_cls and got["subtype"] == want_sub,
          f"got {got['class']}/{got['subtype']} t*={got['t_star']} R_m={got['R_m']}")
# regression against the real anchors' itineraries is covered by T4/S0a.

# ---------------------------------------------------------------- T4
rows = step0()
check("T4 Step-0 gate passes on committed S1 CSVs (all anchors)",
      all(r["pass"] for r in rows),
      "classes: " + ",".join(f"K{r['K']}={r['S0a_class']}" for r in rows))

print()
if FAILURES:
    print(f"FAILING: {FAILURES}")
    raise SystemExit(1)
print("all S2 microscope assertions passed")
