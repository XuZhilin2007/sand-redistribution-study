# -*- coding: utf-8 -*-
"""v04-A second wave: A-W conditional all-time certificate census.

For each (K, alpha) in the probe set:
  1. EXACT steps up to entry t0 = first selection in {l, m}
     (expected: 2 for K>=7 inside the window, 1 for K=6).
  2. EXACT prefix dominance: U(t) >= U(t0) for t < t0.
  3. EXACT square-completion conditions S1'-S4' anchored at t0
     (B-S certificate): S1 A>0, S2 H2>0, S3 H'(1-alpha)<=0,
     S4 H(1-alpha)>=0.
  4. FLOAT residence window check (T=300): selections stay in {l,m}
     after t0; near-ties flagged by margin.
Classification:
  CERTIFIED  = prefix ok AND S1'-S4' ok  (all-time GIVEN residence)
  S-FAIL     = some condition fails (list which)
  PREFIX-FAIL, ENTRY-UNUSUAL recorded.
Validation layer: alpha=1/4 must reproduce the v0.3.0 picture
  (S-FAIL exactly on the finite exception set {9,10,11,16,17,23} for
  CYCLE classes; LOCK classes certified), modulo the anchor-time detail.

Output: cert_census.json + console summary. Exit 0.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

import numpy as np

KS = list(range(6, 204)) + [300, 600, 999, 1000, 1001]
ALPHAS = [F(1, 4), F(247, 1000), F(27, 100), F(29, 100)]
T_FLOAT = 300
FINITE_EXCEPTIONS = {9, 10, 11, 16, 17, 23}


def G(x: F) -> F:
    if x <= 0:
        return F(0)
    if x < F(1, 3):
        return x / 2
    if x <= F(2, 3):
        return 3 * x / 2 - F(1, 3)
    if x <= F(5, 6):
        return 2 * x - F(2, 3)
    return F(1)


def exact_entry(K: int, alpha: F, max_t: int = 12):
    """Run exact steps until first boundary selection; return
    (states, sels, t0)."""
    p0 = [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]
    g = [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]
    m = -(-5 * K // 6)
    states, sels = [p0], []
    t0 = None
    while len(sels) < max_t:
        p = states[-1]
        acc, Dbest, best = F(0), None, 0
        for j, pj in enumerate(p):
            acc += pj
            d = acc - F(j + 1, K)
            if Dbest is None or d > Dbest:
                Dbest, best = d, j
        s = best + 1
        M = alpha * sum(p[:s])
        states.append([((1 - alpha) * pj if j < s else pj) + M * gj
                       for j, (pj, gj) in enumerate(zip(p, g))])
        sels.append(s)
        if s in (m - 1, m):
            # anchor = the state FROM which the first boundary selection
            # is made (E semantics: t0 = first-boundary-selection index)
            t0 = len(sels) - 1
            break
    return states, sels, t0, m


def certificate(K: int, alpha: F, states, t0, m):
    """Exact S1'-S4' at anchor t0. Returns (ok, details, margins)."""
    h = F(1, K)
    b = F(m, K)
    c = 2 * b - b * b
    q = F(6 - K % 6, 3 * K)
    g = [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]
    p0 = states[0]
    U = [K * sum((y - h) ** 2 for y in st) for st in states[:t0 + 1]]
    prefix_ok = all(U[t] >= U[t0] for t in range(t0))
    x0 = states[t0][m - 1]
    gam = lambda x: (c - x) / (1 - q)
    e = [states[t0][i] - gam(x0) * g[i] for i in range(m - 1)]
    G1 = sum(g[:m - 1])
    G2 = sum(v * v for v in g[:m - 1])
    E1 = sum(e)
    GE = sum(g[i] * e[i] for i in range(m - 1))
    E2 = sum(v * v for v in e)
    Sfx = K * sum((p0[i] - h) ** 2 for i in range(m, K))
    gA, gB = c / (1 - q), 1 / (1 - q)
    A = K * gB * gB * G2 + K
    B = -2 * K * gA * gB * G2 + 2 * gB * G1 - 2
    Cc = (K * gA * gA * G2 - 2 * gA * G1 + F(m - 1, K) + F(1, K) + Sfx)
    d = -K * gB * GE
    f_ = K * gA * GE - E1
    E = K * E2
    H2 = E - d * d / A
    H1 = 2 * f_ - B * d / A
    H0 = Cc - U[t0] - B * B / (4 * A)
    rho = 1 - alpha
    cond = {"S1": A > 0, "S2": H2 > 0,
            "S3": 2 * H2 * rho + H1 <= 0,
            "S4": H2 * rho * rho + H1 * rho + H0 >= 0}
    margins = {"S3_margin": float((-(2 * H2 * rho + H1)) / abs(E)),
               "S4_margin": float((H2 * rho * rho + H1 * rho + H0) / abs(E))}
    return prefix_ok, cond, margins


def residence_float(K: int, alpha: float, t0: int, T: int = T_FLOAT):
    j = np.arange(1, K + 1)
    p = (2 / K - (2 * j - 1) / (K * K)).astype(float)
    edges = np.arange(0, K + 1) / K
    Gv = np.where(edges < 1 / 3, edges / 2,
                  np.where(edges <= 2 / 3, 1.5 * edges - 1 / 3,
                           np.where(edges <= 5 / 6, 2 * edges - 2 / 3, 1.0)))
    g = np.diff(Gv)
    jh = j / K
    m = -(-5 * K // 6)
    min_tie = 1.0
    for t in range(T):
        D = np.cumsum(p) - jh
        dmax = D.max()
        s = int(np.argmax(D)) + 1
        # margin between top and the best non-boundary competitor
        mask = np.ones(K, bool)
        mask[m - 2] = False
        mask[m - 1] = False
        if mask.any():
            min_tie = min(min_tie, float(dmax - D[mask].max()))
        if t >= t0 and s not in (m - 1, m):
            return False, min_tie
        M = alpha * float(np.sum(p[:s]))
        p = np.where(np.arange(K) < s, (1 - alpha) * p, p) + M * g
    return True, min_tie


def main() -> int:
    t_start = time.time()
    rows = []
    for alpha in ALPHAS:
        summary = {"certified": 0, "s_fail": 0, "prefix_fail": 0,
                   "entry_unusual": 0, "residence_flag": 0}
        s_fail_ks = []
        for K in KS:
            states, sels, t0, m = exact_entry(K, alpha)
            if t0 is None:
                summary["entry_unusual"] += 1
                rows.append({"K": K, "alpha": str(alpha), "class":
                             "ENTRY-UNUSUAL"})
                continue
            prefix_ok, cond, margins = certificate(
                K, alpha, states, t0, m)
            res_ok, tie = residence_float(K, float(alpha), t0)
            ok = prefix_ok and all(cond.values())
            if not prefix_ok:
                cls = "PREFIX-FAIL"
                summary["prefix_fail"] += 1
            elif not all(cond.values()):
                cls = "S-FAIL:" + ",".join(
                    k for k, v in cond.items() if not v)
                summary["s_fail"] += 1
                s_fail_ks.append(K)
            else:
                cls = "CERTIFIED" + ("" if res_ok else "|RES-FLAG")
                summary["certified"] += 1
                if not res_ok:
                    summary["residence_flag"] += 1
            rows.append({"K": K, "alpha": str(alpha), "t0": t0,
                         "class": cls, "margins": margins,
                         "residence_float": res_ok,
                         "tie_margin": tie})
        print(f"alpha={alpha} ({float(alpha):.4f}) over {len(KS)} K: "
              f"{summary}")
        print(f"  S-FAIL Ks: {sorted(s_fail_ks)}"
              + ("  (K=6 expected: registered t*=1 special case)"
                 if 6 in s_fail_ks else ""))
    # acceptance gate: alpha=1/4 must reproduce the v0.3.0 structure
    sys.path.insert(0, str(Path(__file__).parent))
    from cert_census_validate import validate_baseline
    baseline_ok = validate_baseline(rows)
    print(f"baseline validation (alpha=1/4 vs v0.3.0 structure): "
          f"{'PASS' if baseline_ok else 'FAIL'}")
    dest = Path(__file__).with_name("cert_census.json")
    dest.write_text(json.dumps(rows), encoding="utf-8")
    print(f"\n{len(rows)} pairs in {time.time() - t_start:.1f}s -> {dest}")
    return 0 if baseline_ok else 1


if __name__ == "__main__":
    sys.exit(main())
