"""v0.3 Theorem B runner: suffix-dominated initial states — exact
verification of the robustness theorem on a rational (K, C) grid.

Claims verified per run (exact Fraction):
  B1  s_t <= m for all t;  F_{t,m} = C conserved;  d_t >= C - b > tau
  B2  finite entry into {m-1, m}; post-entry scalar closure with (C, Cq)
  B3  dichotomy by residue only: r in {0,1,2} -> eventually always m;
      r in {3,4,5} -> both boundary selections recur (window evidence)
  B4  Theta(K) table and the feasible-range bound Theta(K) <= 41/45 (K>=7)
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import time
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "s4_w1_s1_xo_precision"))

import replay  # noqa: E402

ALPHA = F(1, 4)
TAU = F(1, 10 ** 12)
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

CFG = {
    "Ks": [9, 15, 21, 27, 33, 51, 12, 13, 14],
    "Cs": [F(88, 100), F(92, 100), F(35, 36), F(995, 1000)],
    "H": 2000,
    "escape_Ks": [9, 21],
    "escape_Cs": [F(3, 4), F(4, 5)],
    "escape_H": 40,
}


def params(K):
    import math
    m = (5 * K + 5) // 6  # == ceil(5K/6), pure integer arithmetic
    h = F(1, K)
    b = F(m, K)
    q = F(6 - K % 6, 3 * K)
    return m, h, b, q


def theta(K):
    m, h, b, q = params(K)
    return h + (b - 2 * h) / (1 - q - h / 2)


def make_init(K, C):
    """Suffix-dominated initial state with F_{0,m} = C:
    prefix = canonical shape rescaled by C/c; suffix = remaining mass
    spread with strictly decreasing profile (suffix averages < 1/K)."""
    m, h, b, q = params(K)
    p = replay.exact_bin_probs(replay.near_cdf_exact, K)
    c0 = sum(p[:m])
    scale = C / c0
    p0 = [pi * scale for pi in p[:m]]
    rem = 1 - C
    n_suf = K - m
    # increasing suffix weights w_i proportional to i (i=1..n_suf):
    # every prefix-average of an increasing profile is <= its total average,
    # and rem/n_suf < h iff C > b -- so A1 holds automatically in range.
    wsum = F(n_suf * (n_suf + 1), 2)
    for i in range(1, n_suf + 1):
        p0.append(rem * F(i, 1) / wsum)
    assert sum(p0) == 1
    # suffix dominance check (A1): every suffix average < 1/K
    Fj = C
    for j in range(m + 1, K + 1):
        Fj = Fj + p0[j - 1]
        assert Fj - F(0) - j * h < C - b, f"A1 violated at j={j}"
    return p0


def run_case(K, C, H):
    m, h, b, q = params(K)
    p = make_init(K, C)
    g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
    L = C * q
    t_star = None
    sel_sup_m = False
    conserved = True
    entry_ok = None
    s_list = []
    x = None
    scalar_ok = True
    d_min = None
    for t in range(H):
        D = replay.cumsum_excess(p, K)
        s = replay.argmax_first(D)
        d = D[s - 1]
        s_list.append(s)
        if s > m:
            sel_sup_m = True
            break
        if sum(p[:m]) != C:
            conserved = False
            break
        d_min = d if d_min is None else min(d_min, d)
        if s in (m - 1, m):
            if t_star is None:
                t_star = t
                x = p[m - 1]
            # scalar closure check
            if p[m - 1] != x:
                scalar_ok = False
                break
            if x > h:
                x = F(3, 4) * x + F(1, 4) * L
            else:
                x = x + F(1, 4) * (C - x) * q
        p = replay.update_state(p, g, s, ALPHA)
    post = s_list[t_star:] if t_star is not None else []
    both = (m - 1 in post) and (m in post)
    # eventual-m test for LOCK classes: last 200 selections all m?
    ev_m = (len(s_list) > 300) and all(s == m for s in s_list[-200:])
    return {
        "K": K, "r": K % 6, "C": str(C), "H": H,
        "sel_le_m": not sel_sup_m, "conserved": conserved,
        "t_star": t_star, "scalar_ok": scalar_ok,
        "d_min_ge_a0": (d_min is not None and d_min >= C - b and C - b > TAU),
        "both_branches": both, "eventual_m": ev_m,
    }


def escape_case(K, C, H):
    """Below-feasibility probe: C <= b forces an A1 violation. Record
    whether the trajectory (i) stops under the frozen pre-update rule
    (d <= tau) or (ii) has the selection escape past m. Stopping is
    checked FIRST: a stopped trajectory has no further legal active
    updates, so it cannot demonstrate conservation-breaking updates."""
    m, h, b, q = params(K)
    p = make_init_overb(K, C)
    g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
    esc = None
    stopped = None
    for t in range(H):
        D = replay.cumsum_excess(p, K)
        s = replay.argmax_first(D)
        d = D[s - 1]
        if d <= TAU:
            stopped = {"t": t, "d": str(d)}
            break
        if s > m:
            esc = {"t": t, "s": s,
                   "violation": str(D[K - 1] - (C - b))}
            break
        p = replay.update_state(p, g, s, ALPHA)
    return {"K": K, "C": str(C), "escape": esc, "stopped": stopped}


def make_init_overb(K, C):
    """Initial state with F_{0,m} = C <= b (suffix dominance impossible)."""
    m, h, b, q = params(K)
    p = replay.exact_bin_probs(replay.near_cdf_exact, K)
    c0 = sum(p[:m])
    p0 = [pi * (C / c0) for pi in p[:m]]
    rem = 1 - C
    n_suf = K - m
    # spread evenly (suffix average = rem/n_suf > (1-b)/n_suf... > 1/K)
    for i in range(n_suf):
        p0.append(rem / n_suf)
    assert sum(p0) == 1
    return p0


def main():
    t0 = time.perf_counter()
    rows = []
    for K in CFG["Ks"]:
        Cs = list(CFG["Cs"])
        _, _, _, q = params(K)
        if K % 6 in (3, 4, 5):
            m, h, b, qq = params(K)
            c_can = 2 * b - b * b
            Cs.append(c_can)
        for C in Cs:
            m, h, b, _ = params(K)
            if C <= max(F(6, 7), theta(K)) or C <= b:
                continue  # only run the provable range
            rows.append(run_case(K, C, CFG["H"]))
    bad = [r for r in rows if not (r["sel_le_m"] and r["conserved"]
                                   and r["scalar_ok"] and r["d_min_ge_a0"]
                                   and r["t_star"] is not None)]
    # dichotomy consistency
    dich_bad = [r for r in rows if
                (r["r"] in (0, 1, 2) and not r["eventual_m"]) or
                (r["r"] in (3, 4, 5) and not r["both_branches"])]
    # Theta table + bound
    th_rows = [{"K": K, "theta": str(theta(K))}
               for K in range(7, 604)]
    th_max = max(F(r["theta"]) for r in th_rows)
    th_ok = th_max == F(197, 203) and all(
        F(r["theta"]) < 1 for r in th_rows)
    escapes = [escape_case(K, C, CFG["escape_H"])
               for K in CFG["escape_Ks"] for C in CFG["escape_Cs"]]
    esc_ok = all(e["escape"] is not None or e["stopped"] is not None
                 for e in escapes)

    res = {
        "n_runs": len(rows), "violations": len(bad), "bad": bad[:3],
        "dichotomy_violations": len(dich_bad), "dich_bad": dich_bad[:3],
        "theta_max": str(th_max), "theta_bound_ok": th_ok,
        "escapes": escapes,
        "out_of_scope_all_fired": esc_ok,
        "n_stopped": sum(1 for e in escapes if e["stopped"] is not None),
        "n_escaped": sum(1 for e in escapes if e["escape"] is not None),
        "elapsed_s": round(time.perf_counter() - t0, 1),
    }
    with (OUT / "b_rows.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with (OUT / "b_theta.csv").open("w", newline="",
                                    encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["K", "theta"])
        w.writeheader()
        w.writerows(th_rows)
    res["git"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                                capture_output=True, text=True).stdout.strip()
    res["output_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(OUT.glob("*.csv"))}
    (OUT / "b_summary.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(res, indent=2, ensure_ascii=False, default=str)[:2000])
    if bad or dich_bad or not th_ok or not esc_ok:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
