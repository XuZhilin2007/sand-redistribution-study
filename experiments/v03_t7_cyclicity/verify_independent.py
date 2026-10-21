"""Independent verifier for the r=5 family-cover certificate.

SECOND IMPLEMENTATION, written from the mathematical definitions in
docs/V03_T7_CYCLICITY_PROOF.md sec 6.2(c) only; it imports NOTHING from
run_t7.py or sand_m0.word_cyclicity. It reads the committed certificate
JSON and re-derives, with exact Fraction arithmetic:

  1. sigma_5(K) from the S3 definitions (m, b, c, q, L, sigma);
  2. the envelope constant C5 by DIRECT enumeration: for every lattice
     n in [1, N_DIRECT], check |sigma_5(n) - 35/292| <= C5/n exactly
     (the all-n guarantee is the polynomial certificate in run_t7
     step2d; this independent enumeration cross-checks it on a large
     finite range along a different computational path);
  3. the parameter box, both thresholds, band, limit cycle, island;
  4. the family chain (corner-worst crossing and inclusion slacks);
  5. the BFS covering (uncertainty strip -> union of both branches).

Run:  python experiments/v03_t7_cyclicity/verify_independent.py
Exit 0 + INDEPENDENT VERIFY PASS on stdout iff every field matches.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERT = HERE / "results" / "t7_family_cover_cert.json"
N_DIRECT = 5000
SEQ = [2, 2, 3, 2, 2]
SIG_INF = F(35, 292)


def params5(K: int):
    """S3 boundary parameters, fresh implementation (integer ceil)."""
    m = -(-5 * K // 6)
    h = F(1, K)
    b = F(m, K)
    c = 2 * b - b * b
    q = F(6 - K % 6, 3 * K)
    L = c * q
    return m, h, b, c, q, L


def sigma5(K: int) -> F:
    m, h, b, c, q, L = params5(K)
    return q * (c - L) / (4 * (h - L))


def sigma5_of_n(n: int) -> F:
    return sigma5(6 * n + 5)


def check_envelope(C5: F) -> bool:
    """Direct enumeration: |sigma_5(n) - 35/292| <= C5/n for all
    lattice n in [1, N_DIRECT]."""
    for n in range(1, N_DIRECT + 1):
        if abs(sigma5_of_n(n) - SIG_INF) > C5 / n:
            return False
    return True


def theta2(beta: F, sigma: F) -> F:
    return (1 - sigma * (1 + beta)) / (beta * beta)


def main() -> int:
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    fails = []

    def want(name, got, expected):
        if got != expected:
            fails.append(f"{name}: cert={expected} recomputed={got}")

    # --- 1+2. envelope constant, cross-checked by direct enumeration ---
    C5 = F(cert["C5_per_n"])
    if not check_envelope(C5):
        fails.append("direct envelope enumeration violated")
    want("sigma_inf", cert["sigma_inf"], str(SIG_INF))
    want("sequence", cert["sequence"], SEQ)

    # --- 3. box, thresholds, band, limit cycle, island ---
    K0 = cert["K0"]
    n0 = (K0 - 5) // 6
    beta_lo = 1 - F(1, 12 * K0)
    beta_hi = F(1)
    sig_lo, sig_hi = SIG_INF - C5 / n0, SIG_INF + C5 / n0
    want("beta_lo", cert["beta_lo"], str(beta_lo))
    want("sig_lo", cert["sig_lo"], str(sig_lo))
    want("sig_hi", cert["sig_hi"], str(sig_hi))
    th_lo = theta2(beta_hi, sig_hi)
    th_hi = theta2(beta_lo, sig_lo)
    want("th_lo", cert["th_lo"], str(th_lo))
    want("th_hi", cert["th_hi"], str(th_hi))
    want("strip_width", cert["strip_width"], str(th_hi - th_lo))
    band_hi = F(3, 4) * (beta_hi + sig_hi)
    want("band_hi", cert["band_hi"], str(band_hi))
    # box sanity at a few lattice points (independent recomputation)
    for K in (K0, K0 + 600, K0 + 60000):
        _, _, _, _, _, _ = params5(K)
        beta_K = 1 - F(6 - K % 6, 3 * K) / 4
        s_K = sigma5(K)
        if not (beta_lo <= beta_K <= beta_hi and sig_lo <= s_K <= sig_hi):
            fails.append(f"box violated at K={K}")
    # limit cycle: fixed point of the (2,2,3,2,2) composition, beta=1
    lam, mu = F(1), F(0)
    for mbr in SEQ:
        for _ in range(mbr):
            lam, mu = lam, mu + SIG_INF
        lam, mu = lam * F(3, 4), mu * F(3, 4)
    y0 = mu / (1 - lam)
    ys = [y0]
    for mbr in SEQ:
        y = ys[-1]
        for _ in range(mbr):
            y = y + SIG_INF
        ys.append(F(3, 4) * y)
    if ys[5] != ys[0]:
        fails.append("limit cycle does not close")
    margins = []
    for i, mbr in enumerate(SEQ):
        margins.append(1 - (ys[i] + (mbr - 1) * SIG_INF))
        margins.append((ys[i] + mbr * SIG_INF) - 1)
    dmin = min(abs(v) for v in margins)
    want("dmin_limit", cert["dmin_limit"], str(dmin))
    num, den = cert["w_frac"].split("/")
    w = F(int(num), int(den)) * dmin
    island = [(y - w, y + w) for y in ys[:5]]
    want("island", [list(x) for x in cert["island"]],
         [[str(a), str(b)] for a, b in island])

    # --- 4. family chain with corner-worst slacks ---
    min_cross = None
    min_incl = None
    for j in range(5):
        a, bb = ys[j] - w, ys[j] + w
        for s in range(1, SEQ[j] + 1):
            a = beta_lo * a + sig_lo
            bb = beta_hi * bb + sig_hi
            slack = (1 - bb) if s < SEQ[j] else (a - 1)
            if slack <= 0:
                fails.append(f"family chain crossing fails at phase {j}")
            min_cross = slack if min_cross is None else min(min_cross, slack)
        a, bb = F(3, 4) * a, F(3, 4) * bb
        lo_n, hi_n = ys[(j + 1) % 5] - w, ys[(j + 1) % 5] + w
        for slack in (a - lo_n, hi_n - bb):
            if slack <= 0:
                fails.append(f"family chain inclusion fails at phase {j}")
            min_incl = slack if min_incl is None else min(min_incl, slack)
    want("min_crossing_margin", cert["min_crossing_margin"], str(min_cross))
    want("min_inclusion_margin", cert["min_inclusion_margin"],
         str(min_incl))

    # --- 5. BFS covering (strip -> union of both branch images) ---
    def img(u, v, mbr):
        s_lo = sum(beta_lo ** i for i in range(mbr))
        s_hi = sum(beta_hi ** i for i in range(mbr))
        return (F(3, 4) * (beta_lo ** mbr * u + sig_lo * s_lo),
                F(3, 4) * (beta_hi ** mbr * v + sig_hi * s_hi))

    live = [(F(3, 4), band_hi, 0)]
    n_pop = 0
    max_depth = 0
    while live:
        lo, hi, d = live.pop()
        n_pop += 1
        if d > 200:
            fails.append("independent BFS depth exceeded")
            break
        max_depth = max(max_depth, d)
        if any(l <= lo and hi <= u for l, u in island):
            continue
        pieces = []
        if hi <= th_lo:
            pieces = [(lo, hi, 3)]
        elif lo >= th_hi:
            pieces = [(lo, hi, 2)]
        else:
            if lo < th_lo:
                pieces.append((lo, th_lo, 3))
            if max(lo, th_lo) < min(hi, th_hi):
                pieces.append((max(lo, th_lo), min(hi, th_hi), 2))
                pieces.append((max(lo, th_lo), min(hi, th_hi), 3))
            if hi > th_hi:
                pieces.append((th_hi, hi, 2))
        for (u, v, mbr) in pieces:
            a2, b2 = img(u, v, mbr)
            live.append((a2, b2, d + 1))
    want("n_pop", cert["n_pop"], n_pop)
    want("max_depth", cert["max_depth"], max_depth)

    if fails:
        print("INDEPENDENT VERIFY FAIL:")
        for f in fails:
            print("  -", f)
        return 1
    print("INDEPENDENT VERIFY PASS "
          f"(K0={K0}, depth={max_depth}, pop={n_pop}, "
          f"min_incl={min_incl}, min_cross={min_cross})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
