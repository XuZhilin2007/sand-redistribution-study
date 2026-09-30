"""M1C.1 fixed-uniform kernel ablation tests (docs/M1C1_FIXED_UNIFORM_KERNEL_ABLATION.md).

Run:  python tests/test_m1c1_uniform_kernel.py

Variant U = canonical M1A selection with the canonical initial state (q_near)
and a fixed, state-unaware, deficit-unaware but TARGET-MATCHED respray law
q_uniform[i] = 1/K. It walks the exact same canonical fixed-throw code path
as Variant A (`redistribution="throw"`); only the respray probabilities
differ (`redistribution_dist=DISTRIBUTIONS["uniform"]`).

Verifies:
  1. the uniform law is exact on the canonical bin grid (1/K each, sum 1);
  2. selection invariance: identical decision rows and state at t=0 for A
     and U (same state -> same selection operator, same removal);
  3. U's update identity: per-round added mass is exactly M_removed/K in
     every bin, fully redeposited (fixed-throw path with the uniform law);
  4. trajectory invariants for U at the canonical configuration: mass
     conservation, non-negativity, decisions match the canonical argmax
     rule on U's own states;
  5. deterministic reproducibility and post-stop freezing;
  6. already-uniform state (no sweep) and zero-removed-mass behaviour
     (never-active run leaves the state constant);
  7. regression anchors vs committed results: M1C.1 U stopping times and
     stop-state uniformity; M1C A/B summary integrity; M1A.1/M1B.2 exact
     stops; backward compatibility (no new kwarg == canonical behaviour).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.adaptive import cumulative_excess, run_m1a_deterministic
from sand_m0.model import DISTRIBUTIONS

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


K, T_MAX, ALPHA, TOL = 100, 5000, 0.25, 1e-12
K_LIST = [50, 100, 200, 400]
q_near, q_unif = DISTRIBUTIONS["near"], DISTRIBUTIONS["uniform"]

# --- 1. uniform law exactness -------------------------------------------------------
qu = q_unif.bin_probs(K)
check("q_uniform is exactly 1/K per bin on the canonical grid (<= 1e-15)",
      float(np.abs(qu - 1.0 / K).max()) <= 1e-15)
check("q_uniform sums to 1", abs(float(qu.sum()) - 1.0) <= 1e-15)

# --- canonical runs -----------------------------------------------------------------
run_a = run_m1a_deterministic(q_near, K=K, T_max=T_MAX, alpha=ALPHA, stop_tolerance=TOL)
run_u = run_m1a_deterministic(q_near, K=K, T_max=T_MAX, alpha=ALPHA,
                              stop_tolerance=TOL, redistribution="throw",
                              redistribution_dist=q_unif)
states_u = np.asarray(run_u["states"], dtype=float)
stop_u = run_u["stopped_at"]
rows_u = run_u["rows"]

# --- 2. selection invariance (same state -> same decision) --------------------------
ok0 = all(run_a["rows"][0][key] == rows_u[0][key]
          for key in ["active", "j_t", "a_t", "D_max", "D_margin",
                      "prefix_mass", "removed_mass"])
check("A/U: identical decision row at t=0 (same state -> same selection, "
      "same removal)", ok0)
check("A/U: identical initial state (both start from q_near)",
      np.array_equal(np.asarray(run_a["states"], dtype=float)[0], states_u[0]))

ok_dec = True
for t in range(stop_u + 1):
    D = cumulative_excess(states_u[t], K)
    if int(np.argmax(D)) + 1 != rows_u[t]["j_t"]:
        ok_dec = False
check("U: recorded boundaries match the canonical argmax rule on U's own states",
      ok_dec)

# --- 3. U update identity: constant added vector M_removed/K ------------------------
ok_add = ok_redep = True
for t in range(stop_u):
    r = rows_u[t]
    if not r["active"]:
        continue
    p_minus = states_u[t].copy()
    p_minus[: r["j_t"]] *= (1.0 - ALPHA)
    added = states_u[t + 1] - p_minus
    if float(np.abs(added - r["removed_mass"] / K).max()) > 1e-15:
        ok_add = False
    if abs(added.sum() - r["removed_mass"]) > 1e-12:
        ok_redep = False
check("U: per-round added mass is exactly M_removed/K in every bin "
      "(fixed uniform respray, no deficit information)", ok_add)
check("U: removed mass fully redeposited each round", ok_redep)

# --- 4. trajectory invariants -------------------------------------------------------
check("U: mass conservation at every state (|sum - 1| <= 1e-12)",
      float(np.abs(states_u.sum(axis=1) - 1.0).max()) <= 1e-12)
check("U: non-negativity at every state", float(states_u.min()) >= 0.0)

run_u2 = run_m1a_deterministic(q_near, K=K, T_max=T_MAX, alpha=ALPHA,
                               stop_tolerance=TOL, redistribution_dist=q_unif)
check("U: deterministic reproducibility (bitwise identical states)",
      all(np.array_equal(s1, s2)
          for s1, s2 in zip(states_u, np.asarray(run_u2["states"], dtype=float))))
check("U: state frozen after the stop (post-stop states identical)",
      stop_u is not None and stop_u < T_MAX
      and all(np.array_equal(states_u[stop_u], s) for s in states_u[stop_u + 1:]))

# --- 6. already-uniform state and zero removed mass ---------------------------------
uniform_state = np.full(K, 1.0 / K)
D_uniform = cumulative_excess(uniform_state, K)
check("already-uniform state: D_max = 0 -> the canonical rule never sweeps it",
      abs(float(D_uniform.max())) <= 1e-15)
r_never = run_m1a_deterministic(q_near, K=K, T_max=50, alpha=ALPHA,
                                stop_tolerance=10.0, redistribution_dist=q_unif)
check("zero removed mass: with stop_tolerance above any D_max the controller "
      "never sweeps and the state stays constant",
      r_never["stopped_at"] == 0
      and all(np.array_equal(np.asarray(r_never["states"], dtype=float)[0], s)
              for s in np.asarray(r_never["states"], dtype=float)[1:]))

# --- 7. regression anchors vs committed results ------------------------------------
res = Path(__file__).resolve().parents[1] / "experiments/m1c1_uniform_kernel_ablation/results"
summary_csv = res / "variant_summary.csv"
if summary_csv.exists():
    srows = list(csv.DictReader(summary_csv.open(encoding="utf-8")))
    by = {(r["variant"], int(r["K"])): r for r in srows}
    check("anchor: U exact stop at the committed round for all canonical K "
          "(K-independent t=6)",
          all(by[("U", Kk)]["t_exact_stop"] == "6" for Kk in K_LIST))
    ok_ustop = all(float(by[("U", Kk)]["U_density_at_stop"]) <= 6e-3 for Kk in K_LIST)
    check("anchor: U stop-state U_density in the committed K-robust range "
          "(<= 6e-3, NOT machine-uniform)", ok_ustop)
    check("anchor: U D_max / U_density monotone over the active phase "
          "(0 increasing steps)", all(float(by[("U", Kk)]["share_D_max_increases"]) == 0.0
                                      and float(by[("U", Kk)]["share_delta_U_positive"]) == 0.0
                                      for Kk in K_LIST))
    check("anchor: committed A rows unchanged (t_stop 324/699/1445/2876)",
          all(by[("A", Kk)]["t_exact_stop"] == str(ts) for Kk, ts in
              [(50, 324), (100, 699), (200, 1445), (400, 2876)]))
    check("anchor: committed B rows unchanged (t_stop 3 for all K)",
          all(by[("B", Kk)]["t_exact_stop"] == "3" for Kk in K_LIST))
else:
    check("M1C.1 summary CSV present (run the experiment first)", False)

check("backward compatibility: run without redistribution_dist == canonical "
      "A behaviour (bitwise identical states)",
      all(np.array_equal(s1, s2)
          for s1, s2 in zip(np.asarray(run_a["states"], dtype=float),
                            np.asarray(run_m1a_deterministic(
                                q_near, K=K, T_max=1000, alpha=ALPHA, stop_tolerance=TOL
                            )["states"], dtype=float))))

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1C.1 fixed-uniform kernel tests passed")
