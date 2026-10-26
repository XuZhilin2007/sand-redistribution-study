# -*- coding: utf-8 -*-
"""v04-A float phase scan: for each (K, alpha) simulate the canonical XO
model with general alpha (float64), record window argmin, first boundary
entry t*, and whether U(t) < U(2) for some t>=3 (two-round failure).

EXPLORATORY evidence only: near-ties (relative gap < 1e-9) are flagged
EXACT-RECHECK; critical boundaries must be confirmed with Fractions.

Usage: python -B scan_phase.py [--quick]
Writes scan_phase.json + prints per-alpha failure summaries.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

import numpy as np

T_WINDOW = 250
RECHECK_EPS = 1e-9


def G_np(x):
    return np.where(x < 1 / 3, x / 2,
                    np.where(x <= 2 / 3, 1.5 * x - 1 / 3,
                             np.where(x <= 5 / 6, 2 * x - 2 / 3, 1.0)))


def make_model(K: int):
    j = np.arange(1, K + 1)
    p0 = 2 / K - (2 * j - 1) / (K * K)
    edges = np.arange(0, K + 1) / K
    Gv = G_np(edges)
    g = np.diff(Gv)
    h = 1 / K
    m = -(-5 * K // 6)
    return p0.astype(float), g, h, m


def simulate(K: int, alpha: float, T: int = T_WINDOW):
    p0, g, h, m = make_model(K)
    p = p0.copy()
    jh = np.arange(1, K + 1) / K
    U = [K * float(np.sum((p - h) ** 2))]
    sels = []
    for _ in range(T):
        D = np.cumsum(p) - jh
        idx = int(np.argmax(D))          # first max = smallest index
        s = idx + 1                      # 1-based selection
        F_s = float(np.sum(p[:s]))
        M = alpha * F_s
        mask = np.arange(K) < s
        p = np.where(mask, (1 - alpha) * p, p) + M * g
        U.append(K * float(np.sum((p - h) ** 2)))
        sels.append(s)
    U = np.array(U)
    t_star = next((t for t, s in enumerate(sels) if s in (m - 1, m)), None)
    return U, sels, t_star


def classify(K: int, alpha: float):
    U, sels, t_star = simulate(K, alpha)
    U2 = U[2]
    tail = U[3:]
    tmin = int(np.argmin(U))             # earliest window argmin
    viol = bool(np.any(tail < U2)) or bool(U[1] < U2)
    # relative gap of the strongest violation (positive => U dips below U2)
    rel = float((U2 - float(tail.min())) / U2) if len(tail) else 0.0
    rel1 = float((U2 - U[1]) / U2)
    flag = "EXACT-RECHECK" if 0 < max(rel, rel1) < RECHECK_EPS else ""
    return {
        "K": K, "alpha": alpha, "t_star": t_star, "argmin": tmin,
        "viol": viol, "rel_gap_tail": rel, "rel_gap_t1": rel1,
        "flag": flag,
    }


def main() -> int:
    quick = "--quick" in sys.argv
    # alpha grid: 1/4 both sides, coarse-to-fine rationals
    alphas = []
    for i in range(8, 33):               # 0.2 .. 0.8 step 1/40
        alphas.append(i / 40)
    for j in (1, 2, 3, 4):               # +- j/96 around 1/4
        alphas += [0.25 - j / 96, 0.25 + j / 96]
    for j in (1, 2, 4, 6, 8, 12, 16, 20):  # fine: +- j/480
        alphas += [0.25 - j / 480, 0.25 + j / 480]
    alphas = sorted(set(round(a, 6) for a in alphas if 0.05 < a < 0.95))
    ks = list(range(6, 301))
    if not quick:
        ks += list(range(302, 1001, 2))
        ks += [1200, 2400, 4800]
    t0 = time.time()
    rows = []
    for ai, alpha in enumerate(alphas):
        fails = []
        for K in ks:
            r = classify(K, alpha)
            rows.append(r)
            if r["viol"] or r["flag"]:
                fails.append((K, round(r["rel_gap_tail"], 6),
                              round(r["rel_gap_t1"], 6), r["argmin"], r["flag"]))
        fs = " ".join(f"{k}({rg},{r1},t*={a}){fl and '!' or ''}"
                      for k, rg, r1, a, fl in fails[:14])
        print(f"a={alpha:.6f} fails={len(fails):4d}: {fs}", flush=True)
    dest = Path(__file__).with_name("scan_phase.json")
    dest.write_text(json.dumps(rows), encoding="utf-8")
    print(f"\n{len(rows)} pairs in {time.time()-t0:.1f}s -> {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
