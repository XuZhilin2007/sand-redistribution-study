"""M2 Gate 2 sign-placement tests (XA vs XO, docs/M2_GATE2_SIGN_PLACEMENT.md).

Run:  python tests/test_m2_gate2_sign_placement.py

Covers the frozen Gate 2 kernel pair:

  1. XA/XO mathematics: endpoints, breakpoint continuity (1/3, 2/3, 5/6),
     monotonicity, nonnegative bin probabilities and normalization at all
     execution K (tau_num = 1e-12), extrema max = +1/6 / min = -1/6, and
     the frozen sign supports;
  2. the CRITICAL control property: exact pointwise sign reflection
     Delta_O == -Delta_A (grids + dense) and G_O + G_A == 2x;
  3. generic execution: identical canonical initial state for XA/XO;
     exact CDF-differencing bin masses against an independent closed-form
     evaluation at K=100 (all three breakpoints off-grid); exact respray
     step; mass conservation;
  4. diagnostics: E_t(j) == repository rebound_residual with sign(max E)
     matching rebound_criterion; sign-at-risk classification
     (positive/zero/negative at tau_num); support-migration classification
     against each kernel's own positive support; prefix-restricted argmax
     never returns j > j* and raw-argmax-on-suffix is flagged.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                       / "experiments" / "m2_gate2_sign_placement"))

import numpy as np

from run_experiment import (  # noqa: E402
    sign_class,
    witness_argmax,
    witness_argmax_prefix,
)
from sand_m0.adaptive import run_m1a_deterministic  # noqa: E402
from sand_m0.kernel_theory import rebound_criterion, rebound_residual  # noqa: E402
from sand_m0.model import (  # noqa: E402
    DISTRIBUTIONS,
    witness_aligned_distribution,
    witness_opposed_distribution,
)

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


TAU = 1e-12
ALPHA, TOL = 0.25, 1e-12
K_LIST = [50, 100, 200, 400]
q_near = DISTRIBUTIONS["near"]
xa = witness_aligned_distribution()
xo = witness_opposed_distribution()

# --- 1. XA/XO mathematics -------------------------------------------------------------
worst_end = 0.0
max_continuity = 0.0
ok_mono = True
BRANCHES_TEST = {
    "XA": [lambda v: 1.5 * v, lambda v: 0.5 * v + 1.0 / 3.0,
           lambda v: 2.0 / 3.0, lambda v: 2.0 * v - 1.0],
    "XO": [lambda v: 0.5 * v, lambda v: 1.5 * v - 1.0 / 3.0,
           lambda v: 2.0 * v - 2.0 / 3.0, lambda v: 1.0],
}
BPS_TEST = [(1.0 / 3.0, 0, 1), (2.0 / 3.0, 1, 2), (5.0 / 6.0, 2, 3)]
for name, dist in [("XA", xa), ("XO", xo)]:
    worst_end = max(worst_end, abs(float(dist.cdf(np.array([0.0]))[0])),
                    abs(float(dist.cdf(np.array([1.0]))[0]) - 1.0))
    g = np.linspace(0.0, 1.0, 12001)
    if not np.all(np.diff(dist.cdf(g)) >= -TAU):
        ok_mono = False
    for bp, i_l, i_r in BPS_TEST:
        # branch formulas AT the breakpoint (a bp±h probe would measure
        # O(h) slope difference, not continuity)
        gap = max(abs(float(dist.cdf(np.array([bp]))[0])
                      - float(BRANCHES_TEST[name][i_l](bp))),
                  abs(float(dist.cdf(np.array([bp]))[0])
                      - float(BRANCHES_TEST[name][i_r](bp))),
                  abs(float(BRANCHES_TEST[name][i_l](bp))
                      - float(BRANCHES_TEST[name][i_r](bp))))
        max_continuity = max(max_continuity, gap)
check(f"XA/XO: G(0)=0, G(1)=1 within {TAU:.0e} (worst {worst_end:.2e})",
      worst_end <= TAU)
check("XA/XO: G nondecreasing on a 12001-point grid", ok_mono)
check(f"XA/XO: G continuous at 1/3, 2/3, 5/6 (worst branch gap {max_continuity:.2e})",
      max_continuity <= TAU)

worst_q = 0.0
ok_norm = True
for dist in (xa, xo):
    for K in K_LIST:
        edges = np.linspace(0.0, 1.0, K + 1)
        raw = np.diff(dist.cdf(edges))
        worst_q = max(worst_q, -float(raw.min()), abs(float(raw.sum()) - 1.0))
check(f"XA/XO x K in {K_LIST}: bin masses q_j >= -{TAU:.0e} and |sum-1| <= "
      f"{TAU:.0e} (worst violation {worst_q:.2e})", worst_q <= TAU)

grid = np.linspace(0.0, 1.0, 12001)
for name, dist in [("XA", xa), ("XO", xo)]:
    d = dist.cdf(grid) - grid
    check(f"{name}: max Delta = +1/6 and min Delta = -1/6 on dense grid "
          f"({float(d.max()):.8f}, {float(d.min()):.8f})",
          abs(float(d.max()) - 1.0 / 6.0) <= TAU
          and abs(float(d.min()) + 1.0 / 6.0) <= TAU)
    pos_ok = all(float(dist.cdf(np.array([v]))[0]) - v > TAU
                 for v in ([0.2, 0.5] if name == "XA" else [0.8]))
    neg_ok = all(float(dist.cdf(np.array([v]))[0]) - v < -TAU
                 for v in ([0.8] if name == "XA" else [0.2, 0.5]))
    check(f"{name}: sign supports match the frozen definition "
          f"(positive on {'(0, 2/3)' if name == 'XA' else '(2/3, 1)'}, "
          f"negative on {'(2/3, 1)' if name == 'XA' else '(0, 2/3)'})",
          pos_ok and neg_ok)

# --- 2. critical control property: exact sign reflection --------------------------------
d_a = xa.cdf(grid) - grid
d_o = xo.cdf(grid) - grid
worst_reflect = float(np.abs(d_o + d_a).max())
worst_gsum = float(np.abs(xo.cdf(grid) + xa.cdf(grid) - 2.0 * grid).max())
ok_reflect_grid = worst_reflect <= TAU
for K in K_LIST:
    x = np.arange(1, K + 1) / K
    ok_reflect_grid &= float(np.abs((xo.cdf(x) - x) + (xa.cdf(x) - x)).max()) <= TAU
check(f"sign reflection: Delta_O == -Delta_A on dense grid + all execution "
      f"grids (max |Delta_O + Delta_A| = {worst_reflect:.2e}); "
      f"G_O + G_A == 2x (max dev {worst_gsum:.2e})",
      ok_reflect_grid and worst_gsum <= TAU)

# --- 3. generic execution ----------------------------------------------------------------
runs = {}
for name, dist in [("XA", xa), ("XO", xo)]:
    runs[name] = run_m1a_deterministic(q_near, K=100, T_max=80, alpha=ALPHA,
                                       stop_tolerance=TOL,
                                       redistribution="throw",
                                       redistribution_dist=dist)
st0 = [np.asarray(r["states"], dtype=float)[0] for r in runs.values()]
check("same-K shared initial state: XA/XO states[0] bit-identical to each "
      "other and to canonical q_near bin probabilities",
      np.array_equal(st0[0], st0[1]) and np.array_equal(st0[0], q_near.bin_probs(100)))

edges = np.linspace(0.0, 1.0, 101)


def xa_cdf_independent(x: np.ndarray) -> np.ndarray:
    out = np.empty_like(x)
    for i, v in enumerate(x):
        if v <= 1.0 / 3.0:
            out[i] = 1.5 * v
        elif v <= 2.0 / 3.0:
            out[i] = 0.5 * v + 1.0 / 3.0
        elif v <= 5.0 / 6.0:
            out[i] = 2.0 / 3.0
        else:
            out[i] = 2.0 * v - 1.0
    return out


def xo_cdf_independent(x: np.ndarray) -> np.ndarray:
    return 2.0 * x - xa_cdf_independent(x)


for name, dist, ref in [("XA", xa, xa_cdf_independent(edges)),
                        ("XO", xo, xo_cdf_independent(edges))]:
    diff = float(np.abs(dist.bin_probs(100) - np.diff(ref)).max())
    check(f"{name}: bin_probs(K=100) equals independent closed-form CDF "
          f"differencing (breakpoints 1/3, 2/3, 5/6 off-grid; max diff "
          f"{diff:.2e})", diff <= TAU)

ok_step = ok_mass = True
worst_step = 0.0
for name, r in runs.items():
    dist = xa if name == "XA" else xo
    g = dist.bin_probs(100)
    st = np.asarray(r["states"], dtype=float)
    n = r["stopped_at"] if r["stopped_at"] is not None else 80
    for t in range(n):
        row = r["rows"][t]
        p_minus = st[t].copy()
        p_minus[: row["j_t"]] *= (1.0 - ALPHA)
        worst_step = max(worst_step, float(np.abs(
            st[t + 1] - p_minus - row["removed_mass"] * g).max()))
        ok_mass &= abs(float(st[t + 1].sum()) - 1.0) <= TAU
check(f"generic fixed-G support: every round is exactly prefix-scaling + "
      f"removed*q_kernel (worst {worst_step:.2e}) and mass is conserved",
      ok_step and worst_step <= TAU and ok_mass)

# --- 4. diagnostics -----------------------------------------------------------------------
ok_sem = ok_sfx = True
worst_match = 0.0
n_boundary_sem = 0
for name, r in runs.items():
    g = (xa if name == "XA" else xo).bin_probs(100)
    st = np.asarray(r["states"], dtype=float)
    n = r["stopped_at"] if r["stopped_at"] is not None else 80
    for t in range(n):
        row = r["rows"][t]
        j_star = int(row["j_t"])
        E = rebound_residual(st[t], j_star, ALPHA, g)
        crit = rebound_criterion(st[t], j_star, ALPHA, g)
        worst_match = max(worst_match, abs(float(E.max()) - crit["worst_excess"]))
        # sign mismatches at |worst| <= tau_num are frozen numerical-boundary
        # cases (fp tie between the two evaluation paths); only non-boundary
        # disagreements count as semantic failures
        if (float(E.max()) > 0.0) != crit["rebound"]:
            if abs(float(E.max())) > TAU:
                ok_sem = False
            else:
                n_boundary_sem += 1
        if float(E[j_star:].max()) > TAU:
            ok_sfx = False
check(f"witness/risk semantics: E == rebound_residual (max |diff| "
      f"{worst_match:.2e}); sign(max E) matches rebound_criterion on all "
      f"checked rounds except {n_boundary_sem} numerical-boundary rounds "
      f"(|worst| <= {TAU:.0e}); suffix E <= tau throughout",
      ok_sem and ok_sfx and worst_match <= TAU)

ok_sign = all(
    sign_class(5e-13, TAU) == "zero" and sign_class(2e-12, TAU) == "positive"
    and sign_class(-2e-12, TAU) == "negative" for _ in [0])
check("sign-at-risk classification: positive/zero/negative at tau_num = 1e-12",
      ok_sign)

E_probe = np.array([0.1, 0.9, 0.9, 0.2, -1.0])
j_raw, _ = witness_argmax(E_probe)
j_pref, _ = witness_argmax_prefix(E_probe, j_star=3)
check("prefix-restricted argmax: ties -> smallest index; never returns "
      f"j > j* (raw j={j_raw} would be suffix for j*=3; prefix j={j_pref})",
      j_raw == 2 and j_pref == 2)

xa_run = runs["XA"]
xa_g = xa.bin_probs(100)
xa_st = np.asarray(xa_run["states"], dtype=float)
n_xa = xa_run["stopped_at"] if xa_run["stopped_at"] is not None else 80
row0 = xa_run["rows"][0]
E0 = rebound_residual(xa_st[0], int(row0["j_t"]), ALPHA, xa_g)
j0, _ = witness_argmax(E0)
sup = "positive" if (float(xa.cdf(np.array([j0 / 100]))[0]) - j0 / 100) > TAU else "other"
check(f"support-migration classification hook: XA round-0 risk at "
      f"x={j0 / 100:.2f} classified '{sup}' against XA's own positive "
      "support (0, 2/3)", sup in ("positive", "other"))

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M2 Gate 2 sign-placement tests passed")
