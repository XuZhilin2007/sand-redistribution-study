# -*- coding: utf-8 -*-
"""v04-E: metric comparison on FIXED trajectories (canonical, alpha=1/4).

For K in {9, 12, 65} compute the exact trajectory once, then evaluate the
optimal stopping time of SEVERAL metrics on the SAME states:
  U      = K * sum (p-1/K)^2            (density quadratic, baseline)
  TV     = sum |p - 1/K|                (total variation x2)
  Dmax   = max_j D_{t,j}                (max prefix discrepancy)
  W1bin  = sum_j |F_{t,j} - j/K|        (Wasserstein-1 on the grid)
  Urange = K * (max p - min p)
Report per metric: argmin in window, value at t=2 relative to min
(improvement retained by stopping at 2), and the common near-optimal
time set within 5% of each metric's own minimum.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path

T = 120


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


def trajectory(K: int):
    p0 = [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]
    g = [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]
    alpha = F(1, 4)
    states = [p0]
    for _ in range(T):
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
    return states


def metrics(states, K):
    h = F(1, K)
    U, TV, Dmax, W1, Ur = [], [], [], [], []
    for p in states:
        U.append(K * sum((y - h) ** 2 for y in p))
        TV.append(sum(abs(y - h) for y in p))
        acc = F(0)
        dm, w = F(0), F(0)
        for j, pj in enumerate(p):
            acc += pj
            d = acc - F(j + 1, K)
            if d > dm:
                dm = d
            w += abs(d)
        Dmax.append(dm)
        W1.append(w)
        Ur.append(K * (max(p) - min(p)))
    return {"U": U, "TV": TV, "Dmax": Dmax, "W1": W1, "Urange": Ur}


def f(x):
    return float(x)


def main() -> int:
    out = {}
    for K in (9, 12, 65):
        states = trajectory(K)
        M = metrics(states, K)
        rec = {}
        print(f"== K={K} (same trajectory, T={T}) ==")
        for name, vals in M.items():
            vmin = min(vals)
            t_arg = min(range(len(vals)), key=lambda t: vals[t])
            # relative excess of stopping at 2 vs this metric's optimum
            excess2 = (vals[2] - vmin) / vmin if vmin > 0 else None
            # near-optimal set: within 5% of min (absolute, metric own scale)
            near = [t for t in range(len(vals))
                    if vals[t] <= vmin * F(105, 100)]
            # is the metric ever below its t=2 value after t=2?
            later_better = any(vals[t] < vals[2] for t in range(3, len(vals)))
            rec[name] = {"argmin": t_arg,
                         "excess2": None if excess2 is None else f(excess2),
                         "near5pct": [near[0], near[-1], len(near)],
                         "later_better_than_2": later_better,
                         "val2_over_min": f(vals[2] / vmin)}
            print(f"  {name:>6}: argmin={t_arg:>3} "
                  f"val(2)/min={f(vals[2] / vmin):8.4f} "
                  f"later_better={later_better} "
                  f"near5%=[{near[0]},{near[-1]}]x{len(near)}")
        out[K] = rec
    dest = Path(__file__).with_name("metrics_probe.json")
    dest.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
