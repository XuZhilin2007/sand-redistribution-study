# -*- coding: utf-8 -*-
"""v04-A R*-census, committed replacement of the earlier inline probe.

Scope (explicit, review R7): this is a FIXED-alpha K-scan.  It does NOT
establish the alpha->0 limit at fixed K by itself; that requires the
cross-alpha comparison included here (R(K, alpha) for three alphas at
the same K) plus, ultimately, a proof.  Natural stopping is honoured
(canonical model: d >= a > tau always, so it never fires here).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

T = 1200
KS = list(range(12, 49)) + [53, 65, 101, 125]
ALPHAS = [1 / 48, 1 / 96, 1 / 192]


def G_np(x):
    return np.where(x < 1 / 3, x / 2,
                    np.where(x <= 2 / 3, 1.5 * x - 1 / 3,
                             np.where(x <= 5 / 6, 2 * x - 2 / 3, 1.0)))


def sim(K: int, alpha: float, T: int = T):
    j = np.arange(1, K + 1)
    p = (2 / K - (2 * j - 1) / (K * K)).astype(float)
    g = np.diff(G_np(np.arange(0, K + 1) / K))
    h = 1 / K
    jh = j / K
    U = [K * float(np.sum((p - h) ** 2))]
    Fs = []
    for _ in range(T):
        D = np.cumsum(p) - jh
        if float(D.max()) <= 1e-12:
            break
        s = int(np.argmax(D)) + 1
        F_s = float(np.sum(p[:s]))
        Fs.append(F_s)
        M = alpha * F_s
        p = np.where(np.arange(K) < s, (1 - alpha) * p, p) + M * g
        U.append(K * float(np.sum((p - h) ** 2)))
    return np.array(U), Fs


def main() -> int:
    out = []
    print("fixed-alpha K-scan; R = cumulative recovered mass at argmin")
    for K in KS:
        row = {"K": K, "r": K % 6}
        line = f"K={K:>3} r={K % 6}"
        for alpha in ALPHAS:
            U, Fs = sim(K, alpha)
            t_arg = int(np.argmin(U))
            R = sum(alpha * Fs[t] for t in range(t_arg))
            row[f"a={alpha}"] = {"argmin": t_arg, "a*t": alpha * t_arg,
                                 "R": R}
            line += f"  | a={alpha}: t*={t_arg:>3} at={alpha * t_arg:.3f} R={R:.4f}"
        out.append(row)
        print(line)
    # cross-alpha stability of R at fixed K (order-separation evidence)
    drift = [max(r[f"a={a}"]["R"] for a in ALPHAS) -
             min(r[f"a={a}"]["R"] for a in ALPHAS) for r in out]
    print(f"\ncross-alpha R spread at fixed K: max={max(drift):.4f} "
          f"mean={sum(drift) / len(drift):.4f} over {len(drift)} K")
    dest = Path(__file__).with_name("rstar_census.json")
    dest.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
