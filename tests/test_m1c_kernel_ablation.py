"""M1C correction-kernel ablation tests (docs/M1C_CORRECTION_KERNEL_ABLATION.md).

Run:  python tests/test_m1c_kernel_ablation.py

Variant B = canonical M1A selection (cumulative-excess argmax, alpha removal,
same stopping rules) with ONLY the redistribution kernel replaced: the
removed mass is placed back in proportion to each bin's positive deficit
against the uniform target u_i = 1/K (`deficit_fill_redistribution`).

Verifies:
  1. kernel unit behaviour: zero-removed identity, single-bin deficit,
     degenerate equality (removed == total deficit fills exactly to target),
     tiny floating-point deficits;
  2. bug guards: zero total deficit or removed > total deficit raise instead
     of silently falling back to the canonical q kernel;
  3. trajectory invariants for Variant B at the canonical configuration:
     mass conservation, non-negativity, deficit-only redistribution (bins at
     or above target receive exactly zero), no overfill beyond floating-point
     saturation, removed mass fully redeposited, and the deficit identity
     total_deficit = post-removal positive excess + removed mass;
  4. shared selection: identical decision row at t=0 for A and B from the
     same state, and B's recorded decisions match the canonical argmax rule
     applied to B's own states;
  5. deterministic reproducibility and post-stop freezing;
  6. regression anchors vs the committed M1C experiment results (when
     present): Variant B exact stopping at t=3 and machine-precision uniform
     post-stop state, Variant A canonical behaviour unchanged by the new
     keyword (default == explicit "throw").
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.adaptive import (
    cumulative_excess,
    deficit_fill_redistribution,
    run_m1a_deterministic,
)
from sand_m0.model import DISTRIBUTIONS

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


K, T_MAX, ALPHA, TOL = 100, 2000, 0.25, 1e-12
q = DISTRIBUTIONS["near"]

# --- 1. kernel unit behaviour -------------------------------------------------------
u10 = np.full(10, 0.1)
out = deficit_fill_redistribution(u10.copy(), 10, 0.0)
check("kernel: zero removed mass returns the state unchanged (no-op)",
      np.array_equal(out, u10))

p = np.full(10, 0.1)
p[0] = 0.15   # positive excess
p[5] = 0.05   # the only deficit bin (deficit 0.05)
out = deficit_fill_redistribution(p.copy(), 10, 0.03)
check("kernel: one-bin deficit receives the whole removed mass, stays <= target",
      abs(out[5] - 0.08) <= 1e-15
      and abs(out.sum() - (p.sum() + 0.03)) <= 1e-15  # kernel adds on top of p_minus
      and out[0] == 0.15 and out[5] <= 0.1 + 1e-15)
out = deficit_fill_redistribution(p.copy(), 10, 0.05)  # degenerate equality case
check("kernel: removed == total deficit fills the bin exactly to target",
      abs(out[5] - 0.1) <= 1e-12
      and abs(out.sum() - (p.sum() + 0.05)) <= 1e-12)

p_tiny = np.full(10, 0.1)
p_tiny[5] -= 1e-15
out = deficit_fill_redistribution(p_tiny.copy(), 10, 5e-16)
check("kernel: tiny floating-point deficit handled (conserved, no overfill)",
      abs(out.sum() - (1.0 - 5e-16)) <= 1e-14 and out[5] <= 0.1 + 1e-15
      and out.min() >= 0.0)

# --- 2. bug guards (no silent fallback to q) ----------------------------------------
def _raises(fn) -> bool:
    try:
        fn()
    except RuntimeError:
        return True
    return False

check("kernel: zero total deficit with removed mass > 0 raises",
      _raises(lambda: deficit_fill_redistribution(u10.copy(), 10, 0.05)))
p_bad = np.full(10, 0.1)
p_bad[5] = 0.05  # deficit 0.05, removed 0.5 >> deficit: identity violation
check("kernel: removed mass > total deficit raises",
      _raises(lambda: deficit_fill_redistribution(p_bad.copy(), 10, 0.5)))

# --- 3+4+5. Variant B trajectory invariants at the canonical configuration ----------
run_a = run_m1a_deterministic(q, K=K, T_max=T_MAX, alpha=ALPHA, stop_tolerance=TOL,
                              redistribution="throw")
run_b = run_m1a_deterministic(q, K=K, T_max=T_MAX, alpha=ALPHA,
                              stop_tolerance=TOL, redistribution="deficit_fill")
states_b = np.asarray(run_b["states"], dtype=float)
stop_b = run_b["stopped_at"]
rows_b = run_b["rows"]

check("Variant B: mass conservation at every state (|sum - 1| <= 1e-12)",
      float(np.abs(states_b.sum(axis=1) - 1.0).max()) <= 1e-12)
check("Variant B: non-negativity at every state",
      float(states_b.min()) >= 0.0)

ok_sel0 = all(
    run_a["rows"][0][key] == rows_b[0][key]
    for key in ["active", "j_t", "a_t", "D_max", "D_margin", "prefix_mass", "removed_mass"]
)
check("Variant B: identical decision row at t=0 (same state -> same selection, "
      "same removal)", ok_sel0)

ok_dec = True
for t in range(stop_b + 1):
    D = cumulative_excess(states_b[t], K)
    if int(np.argmax(D)) + 1 != rows_b[t]["j_t"]:
        ok_dec = False
check("Variant B: recorded boundaries match the canonical argmax rule on B's own "
      "states (shared selection logic)", ok_dec)

ok_redep = ok_defonly = ok_nofill = ok_ident = True
worst_saturation = 0.0
for t in range(stop_b):
    r = rows_b[t]
    if not r["active"]:
        continue
    j_t = r["j_t"]
    removed = r["removed_mass"]
    p_minus = states_b[t].copy()
    p_minus[:j_t] *= (1.0 - ALPHA)
    p_next = states_b[t + 1]
    added = p_next - p_minus
    if abs(added.sum() - removed) > 1e-12:
        ok_redep = False
    deficit = np.maximum(1.0 / K - p_minus, 0.0)
    excess = np.maximum(p_minus - 1.0 / K, 0.0)
    total_deficit = float(deficit.sum())
    # deficit-only: bins at or above target receive exactly zero
    if np.any(added[p_minus >= 1.0 / K - 1e-15] != 0.0):
        ok_defonly = False
    # no overfill: filled bins land at or below target up to fp saturation
    if np.any(p_next > np.maximum(p_minus, 1.0 / K) + 1e-12):
        ok_nofill = False
    worst_saturation = max(worst_saturation, float((p_next - 1.0 / K).max()))
    # deficit identity: total deficit = post-removal positive excess + removed
    if abs(total_deficit - float(excess.sum()) - removed) > 1e-12:
        ok_ident = False
check("Variant B: removed mass fully redeposited each round "
      "(sum of additions == removed mass)", ok_redep)
check("Variant B: deficit-only redistribution (bins at/above target get exactly 0)",
      ok_defonly)
check("Variant B: no overfill (filled bins stay <= target up to fp saturation)",
      ok_nofill, f"worst excess above target = {worst_saturation:.3e}")
check("Variant B: deficit identity total_deficit = E_minus + removed holds "
      "every round (<= 1e-12)", ok_ident)

run_b2 = run_m1a_deterministic(q, K=K, T_max=T_MAX, alpha=ALPHA, stop_tolerance=TOL,
                               redistribution="deficit_fill")
check("Variant B: deterministic reproducibility (bitwise identical states)",
      all(np.array_equal(s1, s2)
          for s1, s2 in zip(states_b, np.asarray(run_b2["states"], dtype=float))))
check("Variant B: state frozen after the stop (post-stop states identical)",
      stop_b is not None and stop_b < T_MAX
      and all(np.array_equal(states_b[stop_b], s) for s in states_b[stop_b + 1:]))

# --- 6. regression anchors vs committed M1C results ---------------------------------
res_dir = Path(__file__).resolve().parents[1] / "experiments/m1c_kernel_ablation/results"
summary_csv = res_dir / "variant_summary.csv"
if summary_csv.exists():
    srows = list(csv.DictReader(summary_csv.open(encoding="utf-8")))
    b_rows = {int(r["K"]): r for r in srows if r["variant"] == "B"}
    a_rows = {int(r["K"]): r for r in srows if r["variant"] == "A"}
    check("anchor: Variant B exact stop at the committed round for all canonical K "
          "(K-independent t=3)",
          all(b_rows[Kk]["t_exact_stop"] == "3" for Kk in [50, 100, 200, 400]))
    ok_uni = all(float(b_rows[Kk]["U_density_at_stop"]) <= 1e-24
                 for Kk in [50, 100, 200, 400])
    check("anchor: Variant B post-stop state is machine-precision uniform "
          "(U_density <= 1e-24)", ok_uni)
    ok_a = all(a_rows[Kk]["t_exact_stop"] == str(ts) for Kk, ts in
               [(50, 324), (100, 699), (200, 1445), (400, 2876)])
    check("anchor: Variant A exact stops reproduce committed M1A.1/M1B.2 values",
          ok_a)
else:
    check("M1C summary CSV present (run the experiment first)", False)

check("backward compatibility: default run_m1a_deterministic == explicit "
      "redistribution='throw' (bitwise identical states)",
      all(np.array_equal(s1, s2)
          for s1, s2 in zip(run_a["states"],
                            np.asarray(run_m1a_deterministic(
                                q, K=K, T_max=500, alpha=ALPHA, stop_tolerance=TOL
                            )["states"], dtype=float))))

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1C kernel-ablation tests passed")
