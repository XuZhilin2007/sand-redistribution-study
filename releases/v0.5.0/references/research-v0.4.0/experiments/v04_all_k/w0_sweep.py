# -*- coding: utf-8 -*-
"""v04 all-K main line, step 2: EXACT finite-complement W0 sweep and
remainder measurement.

For every K in [KLO, KHI] and alpha in {1999/8000, 1/4, 2001/8000}:
  - exact s0 = floor(K/2), exact true s1 = argmax of D1 (exact replay
    of the first update, smallest index) -- no rounding formula assumed;
  - exact U0, U1, U2 (two exact updates);
  - check W0: U0 - U2 > 0 and U1 - U2 > 0 (exact);
  - record R(K, alpha) = K * (U2 - U2_limit(alpha)) using the exact
    continuum polynomials (dual-branch), to measure the EM constant.

Exit 0 iff W0 holds everywhere in the sweep; prints min margins and
max |R|.  JSONL checkpoint per K (resumable).
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from continuum_limit import f2_pieces, pl_sq_int, psub, peval  # noqa

KLO, KHI = 204, 2048
ALPHAS = [F(1999, 8000), F(1, 4), F(2001, 8000)]

_U2LIM = {}
for _up in (False, True):
    _U2LIM[_up] = pl_sq_int(
        [(y0, y1, sl, psub(ic, [F(1)]))
         for y0, y1, sl, ic in f2_pieces(upper=_up)])


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


def w0_case(K: int, alpha: F):
    h = F(1, K)
    p0 = [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]
    g = [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]
    states = [p0]
    for _ in range(2):
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
    U = [K * sum((y - h) ** 2 for y in st) for st in states]
    d02, d12 = U[0] - U[2], U[1] - U[2]
    return d02, d12, U[2]


def main() -> int:
    ckpt = Path(__file__).with_name("w0_sweep.jsonl")
    t0 = time.time()
    # compute missing K, then aggregate over ALL valid records in the
    # file (round-3 fix: resumed rows must count toward the summary and
    # failure total; invalid rows are reported and excluded)
    done = set()
    if ckpt.exists():
        for line in ckpt.read_text(encoding="utf-8").splitlines():
            if line.strip():
                try:
                    done.add(json.loads(line)["K"])
                except json.JSONDecodeError:
                    print("w0_sweep: invalid checkpoint line skipped")
    for K in range(KLO, KHI + 1):
        if K in done:
            continue
        rec = {"K": K}
        for alpha in ALPHAS:
            d02, d12, U2 = w0_case(K, alpha)
            ysel = F(1, 2) + (F(9, 16) if alpha <= F(16, 63)
                              else F(3, 4)) * alpha
            U2lim = peval(_U2LIM[alpha > F(16, 63)], alpha)
            R = K * (U2 - U2lim)
            rec[str(alpha)] = {"d02": str(d02), "d12": str(d12),
                               "ok": d02 > 0 and d12 > 0,
                               "R": str(R)}
        with open(ckpt, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec) + "\n")
        if K % 200 == 0:
            print(f"computed K={K} ({time.time() - t0:.0f}s)", flush=True)
    # aggregate over ALL valid records (fresh + resumed)
    min02 = min12 = None
    maxR = F(0)
    fails = 0
    n_pts = 0
    for line in ckpt.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        for a in rec:
            if a == "K":
                continue
            n_pts += 1
            d02, d12, R = (F(rec[a]["d02"]), F(rec[a]["d12"]),
                           F(rec[a]["R"]))
            if not rec[a]["ok"]:
                fails += 1
                print(f"W0 FAIL (file) K={rec['K']} alpha={a}: "
                      f"d02={d02} d12={d12}")
            if min02 is None or d02 < min02:
                min02 = d02
            if min12 is None or d12 < min12:
                min12 = d12
            if abs(R) > abs(maxR):
                maxR = R
    import hashlib
    src_sha = hashlib.sha256(
        Path(__file__).read_bytes()).hexdigest()[:16]
    print(f"\nsource: w0_sweep.py sha256[:16]={src_sha}; "
          f"KLO={KLO} KHI={KHI} alphas={[str(a) for a in ALPHAS]}")
    print(f"aggregate over {n_pts} (K,alpha) points: "
          f"W0 failures = {fails}")
    print(f"min U0-U2 = {min02} ({float(min02):.8f})")
    print(f"min U1-U2 = {min12} ({float(min12):.8f})")
    print(f"max |K*(U2-U2lim)| = {abs(maxR)} ({float(abs(maxR)):.8f})")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
