# -*- coding: utf-8 -*-
"""v04-C1 long-run leak probe, float64, WITH natural stopping.

Replaces the earlier inline T=3000 probe which illegally continued
past the natural stop (review R2).  The process stops for good when
d = max_j D_j <= 1e-12 (checked BEFORE each update); everything
reported happens on the legal orbit only.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

T_CAP = 5000
TAU = 1e-12


def G_np(x):
    return np.where(x < 1 / 3, x / 2,
                    np.where(x <= 2 / 3, 1.5 * x - 1 / 3,
                             np.where(x <= 5 / 6, 2 * x - 2 / 3, 1.0)))


def sim(K: int, eps: float, T: int = T_CAP):
    j = np.arange(1, K + 1)
    p = (2 / K - (2 * j - 1) / (K * K)).astype(float)
    g = (1 - eps) * np.diff(G_np(np.arange(0, K + 1) / K)) + eps / K
    h = 1 / K
    jh = j / K
    alpha = 0.25
    m = -(-5 * K // 6)
    U = [K * float(np.sum((p - h) ** 2))]
    sels = []
    stopped = None
    d_stop = None
    for t in range(T):
        D = np.cumsum(p) - jh
        dmax = float(D.max())
        if dmax <= TAU:
            stopped, d_stop = t, dmax
            break
        s = int(np.argmax(D)) + 1
        M = alpha * float(np.sum(p[:s]))
        p = np.where(np.arange(K) < s, (1 - alpha) * p, p) + M * g
        U.append(K * float(np.sum((p - h) ** 2)))
        sels.append(s)
    return np.array(U), sels, m, stopped, d_stop


def main() -> int:
    out = []
    for K in (6, 9, 12):
        for eps in (0.1, 0.05, 0.01, 0.001):
            U, sels, m, stopped, d_stop = sim(K, eps)
            t_arg = int(np.argmin(U))
            beyond = sum(1 for s in sels if s > m)
            rec = {"K": K, "eps": eps, "stopped_at": stopped,
                   "d_stop": d_stop, "steps": len(sels),
                   "argmin": t_arg, "U2": float(U[2]),
                   "U_min": float(U.min()), "U_at_stop": float(U[-1]),
                   "beyond_m_before_stop": beyond,
                   "later_better_than_2": bool(U[3:].min() < U[2])
                   if len(U) > 3 else False}
            out.append(rec)
            print(f"K={K:>2} eps={eps:<6} stop@{stopped} d={d_stop:.2e} "
                  f"steps={len(sels):>4} argmin={t_arg} "
                  f"U2={rec['U2']:.5f} Ustop={rec['U_at_stop']:.5f} "
                  f"beyond={beyond} later<U2={rec['later_better_than_2']}")
    dest = Path(__file__).with_name("leak_longrun.json")
    dest.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
