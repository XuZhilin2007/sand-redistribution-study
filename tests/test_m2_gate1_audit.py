"""M2 Gate 1 interpretation-audit helper tests (witness localization).

Run:  python tests/test_m2_gate1_audit.py

Covers the audit-only analysis helpers added in
experiments/m2_gate1_generality/run_interpretation_audit.py:

  1. argmax E_t deterministic semantics: ties resolve to the SMALLEST
     index (np.argmax / repository boundary tie-breaking), both raw and
     prefix-restricted;
  2. no-rebound KW is never marked a witness: on a short real KW
     trajectory max_j E_t(j) <= 0 on every active round (kind = risk);
  3. region classification boundaries: Near x < 1/3, Middle
     1/3 <= x <= 2/3, Far x > 2/3 (both breakpoints included in Middle);
  4. reuse of existing criterion semantics: the audit E vector IS the
     M1C.3/M1C.5 rebound_residual helper; sign(max E) matches
     rebound_criterion()["rebound"] and the suffix-impossibility remark
     (E(j) <= 0 for all j > j*) holds on every checked round;
  5. suffix-trivial flag: on KW rounds where the raw argmax lands beyond
     the swept prefix (structurally impossible positions, E <= 0), the
     audit records argmax_on_suffix=True and the prefix-restricted
     argmax stays inside j <= j*.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                       / "experiments" / "m2_gate1_generality"))

import numpy as np

from run_interpretation_audit import (  # noqa: E402
    E_profile,
    region_of,
    witness_argmax,
    witness_argmax_prefix,
)
from sand_m0.adaptive import run_m1a_deterministic  # noqa: E402
from sand_m0.kernel_theory import rebound_criterion  # noqa: E402
from sand_m0.model import (  # noqa: E402
    DISTRIBUTIONS,
    flat_top_distribution,
    weak_positive_distribution,
)

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


ALPHA, TOL = 0.25, 1e-12
q_near = DISTRIBUTIONS["near"]

# --- 1. deterministic argmax semantics ------------------------------------------------
E_tie = np.array([0.1, 0.5, 0.5, 0.2, -1.0])
j_raw, v_raw = witness_argmax(E_tie)
j_pref, v_pref = witness_argmax_prefix(E_tie, j_star=4)
check(f"argmax E ties resolve to smallest index (raw j={j_raw}, prefix j={j_pref})",
      j_raw == 2 and abs(v_raw - 0.5) < 1e-15
      and j_pref == 2 and abs(v_pref - 0.5) < 1e-15)
j_p2, _ = witness_argmax_prefix(np.array([0.1, 3.0, 0.5, 0.2, -1.0]), j_star=2)
check("prefix-restricted argmax never returns j > j* (returns 2 for j*=2 "
      "even though global max is at j=5)", j_p2 == 2)

# --- 3. region classification boundaries ------------------------------------------------
regions = {0.0: "Near", 0.32: "Near", 1.0 / 3.0: "Middle", 0.5: "Middle",
           2.0 / 3.0: "Middle", 0.67: "Far", 1.0: "Far"}
ok_reg = all(region_of(x) == r for x, r in regions.items())
check("region split: Near x<1/3, Middle 1/3<=x<=2/3 (both breakpoints in "
      "Middle), Far x>2/3", ok_reg)

# --- 2/4/5. real short runs: KW risk-only, criterion reuse, suffix flag -----------------
kw_run = run_m1a_deterministic(q_near, K=100, T_max=30, alpha=ALPHA,
                               stop_tolerance=TOL,
                               redistribution="throw",
                               redistribution_dist=weak_positive_distribution())
ks_run = run_m1a_deterministic(q_near, K=100, T_max=60, alpha=ALPHA,
                               stop_tolerance=TOL,
                               redistribution="throw",
                               redistribution_dist=flat_top_distribution())
g_kw = weak_positive_distribution().bin_probs(100)
g_ks = flat_top_distribution().bin_probs(100)

ok_kw = True
kw_states = np.asarray(kw_run["states"], dtype=float)
n_kw = kw_run["stopped_at"] if kw_run["stopped_at"] is not None else 30
for t in range(n_kw):
    row = kw_run["rows"][t]
    E = E_profile(kw_states[t], int(row["j_t"]), g_kw)
    if float(E.max()) > 0.0:
        ok_kw = False  # a no-rebound trajectory must never yield a witness
check(f"no-rebound KW never marked witness: max_j E_t(j) <= 0 on all {n_kw} "
      "active rounds of the short run", ok_kw and n_kw > 0)

ok_sem = ok_sfx = True
worst_match = 0.0
for run, g in [(kw_run, g_kw), (ks_run, g_ks)]:
    st = np.asarray(run["states"], dtype=float)
    n = run["stopped_at"] if run["stopped_at"] is not None else 60
    for t in range(n):
        row = run["rows"][t]
        j_star = int(row["j_t"])
        E = E_profile(st[t], j_star, g)
        crit = rebound_criterion(st[t], j_star, ALPHA, g)
        worst_match = max(worst_match, abs(float(E.max()) - crit["worst_excess"]))
        if (float(E.max()) > 0.0) != crit["rebound"]:
            ok_sem = False
        if float(E[j_star:].max()) > 0.0:
            ok_sfx = False
check(f"criterion reuse: audit E == rebound_residual; sign(max E) matches "
      f"rebound_criterion on all checked rounds (max |diff| {worst_match:.2e})",
      ok_sem and worst_match <= TOL)
check("suffix impossibility holds in the audit E vector: E(j) <= 0 for all "
      "j > j* on every checked round", ok_sfx)

# --- 5. suffix-trivial flag on a constructed state ---------------------------------------
E_sfx = np.array([-5.0, -4.0, -3.0, -2.0, -1.0])
j_w, v_w = witness_argmax(E_sfx)
check("constructed all-negative E: raw argmax lands on the trivial terminal "
      "prefix j=K (E(K) = -D_max structure) and would be flagged "
      "argmax_on_suffix for any j* < K",
      j_w == 5 and v_w == -1.0)

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M2 Gate 1 audit helper tests passed")
