"""M1A.1 diagnostics tests (docs/M1A_LONG_HORIZON.md).

Run:  python tests/test_m1a_1_diagnostics.py

Verifies:
  1. the discrete cumulative-excess maximum equals the continuous
     one-sided CDF discrepancy sup_x [F(x) - x] (exact, step-CDF argument);
  2. boundary ambiguity diagnostics on hand-built D profiles (adjacent
     flat peak vs far second peak);
  3. regression anchor: the M1A.1 baseline run (same rule/params) reproduces
     the committed M1A canonical trajectory's first 101 rows bit-for-bit;
  4. long-run invariants (mass, non-negativity, stop semantics unchanged);
  5. recurrence probe sanity on a synthetic periodic sequence;
  6. window statistics identity.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from sand_m0.adaptive import run_m1a_deterministic
from sand_m0.diagnostics import (
    boundary_diagnostics,
    recurrence_l1,
    sup_cdf_minus_uniform,
    window_stats,
)
from sand_m0.model import DISTRIBUTIONS

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


K, TOL = 100, 1e-12

# --- 1. discrete D_max == continuous one-sided sup ---------------------------------
for name in ["near", "uniform", "far"]:
    m = DISTRIBUTIONS[name].bin_probs(K)
    D = np.cumsum(m) - np.arange(1, K + 1) / K
    disc = float(D.max())
    cont = sup_cdf_minus_uniform(m, K)
    check(f"q_{name}: discrete D_max == sup_x [F(x)-x] (one-sided KS functional)",
          abs(disc - cont) < 1e-12, f"disc={disc!r} cont={cont!r}")

# also on an arbitrary M1A state (mid-trajectory)
res_short = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=37,
                                  alpha=0.25, stop_tolerance=TOL)
m37 = res_short["states"][23]
D = np.cumsum(m37) - np.arange(1, K + 1) / K
check("arbitrary M1A state: discrete D_max == continuous sup",
      abs(float(D.max()) - sup_cdf_minus_uniform(m37, K)) < 1e-12)

# --- 2. ambiguity diagnostics on hand-built D profiles -------------------------------
def mass_from_D(D: np.ndarray, K: int) -> np.ndarray:
    """Bin masses realizing a prescribed cumulative excess profile D(1..K).

    m_j = D(j) - D(j-1) + 1/K with D(0)=0; total mass = 1 iff D(K)=0.
    """
    D0 = np.concatenate([[0.0], D])
    return np.diff(D0) + 1.0 / K


# wide single peak: D rises to 0.12 at j=19, falls back to 0 at j=100
D1 = np.concatenate([np.linspace(0.0, 0.12, 20)[1:], np.linspace(0.12, 0.0, 82)[1:]])
assert D1.shape == (K,) and abs(D1[-1]) < 1e-12
d = boundary_diagnostics(mass_from_D(D1, K), K)
check("narrow single peak: one component; top-2 adjacent; A_5% compact",
      d["components_5pct"] == 1 and d["spatial_separation"] <= 0.02
      and d["main_width_5pct"] <= 0.10,
      f"comps5={d['components_5pct']} sep={d['spatial_separation']:.3f} "
      f"main={d['main_width_5pct']:.3f}")

# two far-apart near-equal peaks: 0.10 at j=24, 0.099 at j=82
D2 = np.concatenate([
    np.linspace(0.0, 0.10, 25)[1:],                       # rise to peak 1 (j=24)
    np.linspace(0.10, 0.09, 35)[1:],                      # shallow valley
    np.linspace(0.09, 0.099, 25)[1:],                     # rise to peak 2 (j=82)
    np.linspace(0.099, 0.0, 19)[1:],                      # fall to D(K)=0
])
assert abs(D2[-1]) < 1e-12 and D2.shape == (K,)
d = boundary_diagnostics(mass_from_D(D2, K), K)
check("far second peak: A_1% splits into 2 components with a far member",
      d["components_1pct"] == 2 and d["farthest_point_1pct"] > 0.4,
      f"comps1={d['components_1pct']} far={d['farthest_point_1pct']:.3f}")
check("far second peak: top-2-by-value is the adjacent grid point (grid effect)",
      d["spatial_separation"] <= 0.02,
      f"sep={d['spatial_separation']:.3f}")
check("far second peak: A_5% also split in 2 components",
      d["components_5pct"] == 2 and d["farthest_point_5pct"] > 0.4,
      f"comps5={d['components_5pct']} far5={d['farthest_point_5pct']:.3f}")

# adjacent near-tie: D peaks 0.10 at j=24 and 0.099999 at j=25
D3 = np.concatenate([
    np.linspace(0.0, 0.10, 25)[1:],
    np.array([0.099999, 0.0990]),
    np.linspace(0.0990, 0.0, 75)[1:],
])
assert D3.shape == (K,) and abs(D3[-1]) < 1e-12
d = boundary_diagnostics(mass_from_D(D3, K), K)
check("adjacent near-tie: tiny value gap AND separation 1/K (grid effect, not ambiguity)",
      d["value_gap"] < 1e-5 and d["spatial_separation"] == 0.01,
      f"gap={d['value_gap']:.2e} sep={d['spatial_separation']}")
check("all-D-negative state: diagnostics degrade gracefully (no near-optimal set)",
      boundary_diagnostics(np.zeros(K), K)["components_1pct"] == 1
      and boundary_diagnostics(np.zeros(K), K)["span_1pct"] == 0.0)

# --- 3. regression anchor vs committed M1A canonical trajectory ----------------------
csv_path = Path(__file__).resolve().parents[1] / "experiments/m1a_adaptive_boundary/results/m1a_baseline_trajectory.csv"
m1a_rows = list(csv.DictReader(csv_path.open(encoding="utf-8")))
long_run = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=5000,
                                 alpha=0.25, stop_tolerance=TOL)
same = all(
    abs(long_run["rows"][t]["U_L2"] - float(m1a_rows[t]["U_L2"])) == 0.0
    and long_run["rows"][t]["a_t"] == float(m1a_rows[t]["a_t"])
    for t in range(101)
)
check("M1A.1 long run reproduces committed M1A trajectory bit-for-bit (t=0..100)",
      same)
check("M1A.1 long run: still sweeping at t=5000 iff canonical rule says so",
      (long_run["stopped_at"] is None) == all(r["active"] for r in long_run["rows"]))

# --- 4. long-run invariants ----------------------------------------------------------
ok_mass = all(abs(s.sum() - 1.0) < 1e-12 for s in long_run["states"])
ok_nonneg = all(bool(np.all(s >= -1e-15)) for s in long_run["states"])
check("5000-round run: mass conserved", ok_mass)
check("5000-round run: masses non-negative", ok_nonneg)
check("5000-round run: T_max+1 states recorded",
      len(long_run["states"]) == 5001 and len(long_run["rows"]) == 5001)

# --- 5. recurrence probe sanity -------------------------------------------------------
periodic = [np.roll(np.eye(4)[0], k % 4) * 1.0 for k in range(40)]
rec = recurrence_l1(periodic, max_lag=8, second_half=True)
check("period-4 synthetic sequence: lag-4 L1 distance is exactly 0",
      rec[4]["min"] == 0.0 and rec[8]["min"] == 0.0,
      f"lag4 min={rec[4]['min']}")
random_like = [
    0.25 + 0.1 * np.array([np.sin(i), np.sin(1.7 * i), np.sin(2.3 * i), np.sin(3.1 * i)])
    for i in range(40)
]
rec2 = recurrence_l1(random_like, max_lag=8, second_half=True)
check("non-periodic synthetic sequence: no exactly-recurring lag",
      all(v["min"] > 1e-12 for v in rec2.values()))

# --- 6. window stats identity ---------------------------------------------------------
w = window_stats(np.array([1.0, 2.0, 3.0, 4.0]))
check("window_stats matches hand values",
      w["min"] == 1.0 and w["max"] == 4.0 and abs(w["mean"] - 2.5) < 1e-15
      and abs(w["std"] - np.std([1, 2, 3, 4])) < 1e-15)

print()
if FAILURES:
    print(f"{len(FAILURES)} test(s) FAILED: {FAILURES}")
    raise SystemExit(1)
print("all M1A.1 diagnostics tests passed")
