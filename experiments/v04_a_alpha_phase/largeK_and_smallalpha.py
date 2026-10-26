# -*- coding: utf-8 -*-
"""v04-A supplementary float probes:
 (1) large-K spot checks INSIDE the candidate common window
     (alpha_-, alpha_+) = (0.2461548, 0.2916454): any violation (window
     argmin != 2 or dip below U2) would shrink the window.
 (2) alpha -> 0 organisation: argmin scaling, alpha*argmin, cumulative
     recovered mass R_t = sum_{tau<t} alpha*F_{s_tau} at the argmin.

Float64, exploratory; exact confirmation reserved for critical finds.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

T_LARGE = 400
T_SMALL = 1500


def G_np(x):
    return np.where(x < 1 / 3, x / 2,
                    np.where(x <= 2 / 3, 1.5 * x - 1 / 3,
                             np.where(x <= 5 / 6, 2 * x - 2 / 3, 1.0)))


def simulate(K: int, alpha: float, T: int):
    j = np.arange(1, K + 1)
    p = (2 / K - (2 * j - 1) / (K * K)).astype(float)
    edges = np.arange(0, K + 1) / K
    g = np.diff(G_np(edges))
    h = 1 / K
    jh = j / K
    m = -(-5 * K // 6)
    U = [K * float(np.sum((p - h) ** 2))]
    sels, Fs_list = [], []
    for _ in range(T):
        D = np.cumsum(p) - jh
        idx = int(np.argmax(D))
        s = idx + 1
        F_s = float(np.sum(p[:s]))
        Fs_list.append(F_s)
        M = alpha * F_s
        mask = np.arange(K) < s
        p = np.where(mask, (1 - alpha) * p, p) + M * g
        U.append(K * float(np.sum((p - h) ** 2)))
        sels.append(s)
    return np.array(U), sels, Fs_list, m


def main() -> int:
    out = {"largeK": [], "small_alpha": []}
    print("== (1) large-K inside-window spot checks (T=400) ==")
    # residues {0,2,4} via even anchors; {1,3,5} via odd anchors (R7)
    for K in (5000, 7000, 10000, 14000, 20000, 30000,
              5001, 7001, 10001, 14001, 20001, 30001):
        for alpha in (0.2465, 0.25, 0.27, 0.29, 0.2916):
            U, sels, _, m = simulate(K, alpha, T_LARGE)
            U2 = U[2]
            tail_min = U[3:].min()
            t_arg = int(np.argmin(U))
            viol = bool(tail_min < U2) or bool(U[1] < U2)
            t_star = next((t for t, s in enumerate(sels)
                           if s in (m - 1, m)), None)
            rec = {"K": K, "alpha": alpha, "viol": viol, "argmin": t_arg,
                   "rel_gap": float((U2 - tail_min) / U2), "t_star": t_star}
            out["largeK"].append(rec)
            print(f"K={K:>6} a={alpha:<7} viol={viol} argmin={t_arg} "
                  f"t*={t_star} relgap={rec['rel_gap']:+.2e}")
    print("\n== (2) alpha->0 organisation (float, T=1500) ==")
    for K in (9, 12, 65):
        for alpha in (1 / 16, 1 / 24, 1 / 32, 1 / 48, 1 / 64, 1 / 96):
            U, sels, Fs, m = simulate(K, alpha, T_SMALL)
            t_arg = int(np.argmin(U))
            # cumulative recovered mass up to argmin
            R = sum(alpha * Fs[t] for t in range(t_arg))
            rec = {"K": K, "alpha": alpha, "argmin": t_arg,
                   "alpha_times_argmin": alpha * t_arg, "R_at_argmin": R,
                   "t_star": next((t for t, s in enumerate(sels)
                                   if s in (m - 1, m)), None),
                   "U_arg_over_U2": float(U[t_arg] / U[2])}
            out["small_alpha"].append(rec)
            print(f"K={K:>3} a={alpha:.5f} argmin={t_arg:>4} "
                  f"a*t={alpha * t_arg:.3f} R={R:.4f} "
                  f"t*={rec['t_star']} U*/U2={rec['U_arg_over_U2']:.4f}")
    dest = Path(__file__).with_name("largeK_smallalpha.json")
    dest.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
