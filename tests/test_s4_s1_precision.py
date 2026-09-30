"""S4 Wave 1 / S1 high-precision replay tests (docs/S4_S1_XO_PRECISION.md).

Run:  python tests/test_s4_s1_precision.py

Covers the S1 replay engine (experiments/s4_w1_s1_xo_precision/replay.py):

  1. exact rational bin probabilities: telescoping sum == 1 exactly,
     nonnegativity, agreement with the float implementation's bin_probs,
     and the K=12 on-grid breakpoint value (q_4 = G(1/3) - G(1/4) = 1/24);
  2. the invariant anchor m = ceil(5K/6), b = m/K for the tested K;
  3. K=6 XO hand-check: canonical initial D profile (5,8,9,8,5,0)/36,
     s_0 = 3, d_0 = 1/4, and the proved invariant D_{t,m} = b(1-b) =
     5/36 with s_t <= m over an exact replay slice;
  4. XO K=50 round 0 against the committed Gate-2 float record
     (j_t = 25, D_max = 1/4, moved = 3/16);
  5. the exact worst-excess identity max_j E_t(j) = d_{t+1} - d_t
     (consequence of the M1C.2 decomposition) with Fraction equality,
     and exact tie flag == exact d-recurrence;
  6. mp50/mp100 agreement with the exact tier on an XO K=50 slice
     (selection sequences identical; d_t to <= 1e-38 / <= 1e-55);
  7. cross-check of the diagnostics composition against the repository
     float helpers (rebound_residual / target_matched_excess);
  8. KC exact replay at K=50 stops within the horizon (t_stop recorded,
     compared in the round report against the float 324).
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "s4_w1_s1_xo_precision"))

import numpy as np  # noqa: E402

import replay  # noqa: E402
from sand_m0.kernel_theory import rebound_residual, target_matched_excess  # noqa: E402
from sand_m0.model import DISTRIBUTIONS, witness_opposed_distribution  # noqa: E402

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    if ok:
        print(f"[PASS] {name}" + (f" — {detail}" if detail else ""))
    else:
        print(f"[FAIL] {name} — {detail}")
        FAILURES.append(name)


# ---------------------------------------------------------------- 1
q_near_exact = replay.exact_bin_probs(replay.near_cdf_exact, 50)
q_xo_exact = replay.exact_bin_probs(replay.xo_cdf_exact, 50)
check("exact q_near K=50: sum == 1 exactly, all >= 0",
      sum(q_near_exact) == 1 and all(q >= 0 for q in q_near_exact))
check("exact q_XO K=50: sum == 1 exactly, all >= 0",
      sum(q_xo_exact) == 1 and all(q >= 0 for q in q_xo_exact))
q_near_float = DISTRIBUTIONS["near"].bin_probs(50)
q_xo_float = witness_opposed_distribution().bin_probs(50)
diff_near = max(abs(float(a) - b) for a, b in zip(q_near_exact, q_near_float))
diff_xo = max(abs(float(a) - b) for a, b in zip(q_xo_exact, q_xo_float))
check("exact bin probs match float bin_probs (K=50, both laws)",
      diff_near <= 1e-15 and diff_xo <= 1e-15,
      f"max|diff| near={diff_near:.2e} xo={diff_xo:.2e}")
q12 = replay.exact_bin_probs(replay.xo_cdf_exact, 12)
check("XO K=12 on-grid breakpoint: q_4 == G(1/3)-G(1/4) == 1/24",
      q12[3] == Fraction(1, 6) - Fraction(1, 8) == Fraction(1, 24))

# ---------------------------------------------------------------- 2
expected_m = {6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 10, 50: 42, 100: 84,
              200: 167, 400: 334}
ok_m = all(replay.anchor_m_b(K)[0] == m for K, m in expected_m.items())
check("anchor m = ceil(5K/6) for tested K", ok_m)
b50 = Fraction(42, 50)
b400 = Fraction(334, 400)
check("anchor values b(1-b): 84/625 @ K=50,100 and 5511/40000 @ K=200,400",
      b50 * (1 - b50) == Fraction(84, 625)
      and b400 * (1 - b400) == Fraction(5511, 40000) == Fraction(137775, 1000000))

# ---------------------------------------------------------------- 3
r6 = replay.replay("XO", 6, "exact", 50, identity_check=True)
D0 = replay.cumsum_excess(replay.exact_bin_probs(replay.near_cdf_exact, 6), 6)
check("K=6 initial D profile == (5,8,9,8,5,0)/36",
      D0 == [Fraction(n, 36) for n in (5, 8, 9, 8, 5, 0)])
check("K=6 round 0: s_0 = 3, d_0 = 1/4",
      r6["rows"][0]["s"] == 3 and r6["rows"][0]["_d"] == Fraction(1, 4))
inv6 = all(row["invariant_exact_ok"] for row in r6["rows"])
sb6 = all(row["s_bound_ok"] for row in r6["rows"])
check("K=6 exact slice: invariant D_{t,5} = 5/36 and s_t <= 5 every round",
      inv6 and sb6 and r6["exact_flags_ok"] and r6["identity_ok"],
      f"H=50, censored={r6['censored']}")

# ---------------------------------------------------------------- 4
r50 = replay.replay("XO", 50, "exact", 3, identity_check=True)
row0 = r50["rows"][0]
check("XO K=50 round 0 vs committed float record: j_t=25, d=1/4, moved=3/16",
      row0["s"] == 25 and row0["_d"] == Fraction(1, 4)
      and float(row0["moved"]) == 0.1875)
check("XO K=50 round 0 excess e_0 = 1/4 - 84/625 = 289/2500",
      row0["_d"] - r50["anchor"] == Fraction(289, 2500))

# ---------------------------------------------------------------- 5
r50i = replay.replay("XO", 50, "exact", 40, identity_check=True)
rows = r50i["rows"]
ties = sum(1 for r in rows if r["tie_exact"])
recs = sum(1 for a, b in zip(rows, rows[1:]) if b["_d"] == a["_d"])
D_final = replay.cumsum_excess(r50i["final_state"], 50)
d_after = D_final[replay.argmax_first(D_final) - 1]
if d_after == rows[-1]["_d"]:
    recs += 1
check("exact identity worst == d_{t+1}-d_t holds (K=50, 40 rounds)",
      r50i["identity_ok"])
check("exact tie flag (worst == 0) == exact d-recurrence (incl. final state)",
      ties == recs, f"ties={ties}, d-recurrences={recs}")

# ---------------------------------------------------------------- 6
r50_p50 = replay.replay("XO", 50, "mp50", 40)
r50_p100 = replay.replay("XO", 50, "mp100", 40)

import mpmath as mp  # noqa: E402

mp.mp.dps = 100


def _m(s: str):
    return mp.mpf(s)


d_ex = [_m(r["d"]) for r in rows]
d_50 = [_m(r["d"]) for r in r50_p50["rows"]]
d_100 = [_m(r["d"]) for r in r50_p100["rows"]]
e50 = max(abs(a - b) for a, b in zip(d_ex, d_50))
e100 = max(abs(a - b) for a, b in zip(d_ex, d_100))
sel_ok = ([r["s"] for r in rows] == [r["s"] for r in r50_p50["rows"]]
          == [r["s"] for r in r50_p100["rows"]])
check("mp50/mp100 selection sequences identical to exact (K=50, 40 rounds)",
      sel_ok)
check("mp50 d_t matches exact to <= 1e-38 (parsed at 100 dps)",
      e50 <= mp.mpf(10) ** -38, f"max|diff|={mp.nstr(e50, 3)}")
check("mp100 d_t matches exact to <= 1e-55 (parsed at 100 dps)",
      e100 <= mp.mpf(10) ** -55, f"max|diff|={mp.nstr(e100, 3)}")

# ---------------------------------------------------------------- 7
# Pair each round's diagnostics with its OWN round's state: state_t is the
# final_state of the H=t replay (state_0 = no sweeps yet).
r_h = {h: replay.replay("XO", 50, "exact", h) for h in (0, 1, 2, 3)}
r3 = r_h[3]
qx = replay.exact_bin_probs(replay.xo_cdf_exact, 50)
g_xo = witness_opposed_distribution().bin_probs(50)
worst_diffs = []
for t in (0, 1, 2):
    state_t = r_h[t]["final_state"]
    mass_t = np.array([float(x) for x in state_t])
    j_t = r3["rows"][t]["s"]
    E_repo = rebound_residual(mass_t, j_t, 0.25, g_xo)
    dg = replay.diagnostics(state_t, qx, 50, j_t, Fraction(1, 4))
    worst_diffs.append(abs(float(E_repo.max()) - float(dg["worst"])))
check("diagnostics worst matches repository rebound_residual (t=0,1,2)",
      max(worst_diffs) <= 1e-15,
      f"max|diff|={max(worst_diffs):.2e}")
D_tm_repo = target_matched_excess(
    np.array([float(x) for x in r_h[2]["final_state"]]),
    r3["rows"][2]["s"], 0.25)
state2 = r_h[2]["final_state"]
D2 = replay.cumsum_excess(state2, 50)
j2 = r3["rows"][2]["s"]
gamma_exact = D2[j2 - 1] - (Fraction(3, 4) * D2[j2 - 1]
                            - Fraction(1, 4) * Fraction(j2, 50)
                            * (1 - sum(state2[:j2])))
check("gamma composition matches target_matched_excess (t=2, K=50)",
      abs(float(gamma_exact) - float(D2[j2 - 1]) + float(D_tm_repo[j2 - 1]))
      <= 1e-15)

# ---------------------------------------------------------------- 8
kc = replay.replay("KC", 50, "exact", 500)
check("KC exact K=50 stops within H=500 (float reference t_stop=324)",
      kc["stopped_at"] is not None, f"exact t_stop={kc['stopped_at']}")

print()
if FAILURES:
    print(f"FAILING: {FAILURES}")
    raise SystemExit(1)
print("all S1 precision-engine assertions passed")
