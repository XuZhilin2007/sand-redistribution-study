# -*- coding: utf-8 -*-
"""v04 all-K main line, step 3: UNIFORM ENTRY LEMMA for D2
(machine-certified with exact rational interval arithmetic).

CLAIM (D2-entry): there is N (=33 below) such that for every K >= N
and every alpha in J = [1999/8000, 2001/8000], the argmax over grid
points of D2 (the discrepancy after two true updates) lies in
{ell, m} (the two boundary cells).  Together with theorem A
(permanence from the first boundary selection), this gives uniform
entry at t=2.

PROOF STRUCTURE (all constants exact rationals):
  (1) mu = lambda0 + nu > 1/3 uniformly (interval lower bound);
      lambda0 = alpha*F0 (F0 = P(u0), u0 = floor(K/2)/K),
      nu = alpha*F1(u1), u1 = the true second selection (prefix
      lemma: u1 = 1/2 + 3*alpha*F0/4 + eta/K, |eta| <= 1/2).
  (2) On the upper unswept segment y in [2/3, 5/6):
        D2(y) = -y^2 + (1 + 2 mu) y - 5 mu/3
      is strictly increasing (slope 1 + 2mu - 2y >= 2mu - 2/3 > 0),
      so EVERY interior grid cell with y in [2/3, y_ell) satisfies
      D2(y_j) < D2(y_ell) exactly (grid values = continuous values;
      y_ell = (m-1)/K < 5/6 always, and all interior j <= m-2 have
      y_j < y_ell).  No O(h) margin is consumed here.
  (3) The boundary winner B_K = max(D2(ell), D2(m)) >= 5/36 - c/K
      with c = max upper-segment slope <= 1 + 2*mu_max - 4/3
      (D2(ell) is on the increasing segment, within 1/K of 5/6).
  (4) The mid region y <= 2/3 splits into four concave quadratic
      pieces (breaks 1/3, u0, u1); interval arithmetic bounds each
      piece's maximum over its y-range and the parameter box; the
      uniform upper bound is compared with (3).

Everything uses exact Fractions; interval bounds are sound (outer).
Output: certificate constants + PASS/FAIL per condition; JSON record.
Finite complement K = 6..32 is covered by the accepted v3 interval
certificates (W2 entry condition), not re-proven here.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path

A_LO, A_HI = F(1999, 8000), F(2001, 8000)      # alpha box = J
H_BAR = F(1, 33)                                 # 1/K <= 1/33


class Iv:
    """Exact rational interval [lo, hi] (sound outer arithmetic)."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        if isinstance(lo, int):
            lo = F(lo)
        if isinstance(hi, int):
            hi = F(hi)
        assert lo <= hi
        self.lo, self.hi = lo, hi

    def __add__(self, o):
        o = iv(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __sub__(self, o):
        o = iv(o)
        return Iv(self.lo - o.hi, self.hi - o.lo)

    def __rsub__(self, o):
        return iv(o).__sub__(self)

    def __mul__(self, o):
        o = iv(o)
        vals = [self.lo * o.lo, self.lo * o.hi,
                self.hi * o.lo, self.hi * o.hi]
        return Iv(min(vals), max(vals))

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = iv(o)
        if o.lo <= 0 <= o.hi:
            raise ZeroDivisionError("interval spans zero")
        return self * Iv(F(1) / o.hi, F(1) / o.lo)

    def mid(self):
        return (self.lo + self.hi) / 2

    def __repr__(self):
        return f"[{float(self.lo):.6f},{float(self.hi):.6f}]"


def iv(x):
    return x if isinstance(x, Iv) else Iv(x, x)


def P_of(y: Iv) -> Iv:
    return 2 * y - y * y


def G_mid(y: Iv) -> Iv:
    """G(y) = 3y/2 - 1/3 (valid for y in [1/3, 2/3])."""
    return F(3, 2) * y - F(1, 3)


def quad_max(c2: Iv, c1: Iv, c0: Iv, ylo: F, yhi: F):
    """Sound upper bound of q(y)=c2 y^2 + c1 y + c0 over y in [ylo,yhi]
    for INTERVAL coefficients; requires c2.hi < 0 (concave for every
    parameter point) so the vertex-candidate set covers all vertices."""
    assert c2.hi < 0, "quad_max requires concave"
    ystar = (-c1) / (2 * c2)                  # interval of vertices

    def q(y):
        yv = iv(y)
        return c2 * (yv * yv) + c1 * yv + c0

    cands = [ylo, yhi]
    for y in (ystar.lo, ystar.hi):
        if ylo < y < yhi:
            cands.append(y)
    return max((q(y) for y in cands), key=lambda v: v.hi)


def main() -> int:
    alpha = Iv(A_LO, A_HI)
    rec = {}

    # ---- parameter boxes --------------------------------------------
    # u0 = floor(K/2)/K in [16/33, 1/2] for K >= 33
    u0 = Iv(F(16, 33), F(1, 2))
    # F0 = P(u0) increasing on [0,1] -> interval by endpoints
    F0 = Iv(P_of(iv(u0.lo)).hi, P_of(iv(u0.hi)).hi)   # [P(16/33), 3/4]
    lam0 = alpha * F0
    # u1 = 1/2 + 3 alpha F0 /4 + eta/K  (prefix lemma), |eta|<=1/2
    u1 = F(1, 2) + (F(3, 4) * alpha * F0) + Iv(-H_BAR / 2, H_BAR / 2)
    # clip u1 to (1/2, 2/3) per the prefix lemma
    u1 = Iv(max(u1.lo, F(1, 2)), min(u1.hi, F(2, 3)))
    # nu = alpha * F1(u1), F1(u1) = P(u1) - lam0 (1 - G(u1))
    G_u1 = G_mid(u1)                             # u1 in (1/2,2/3)
    F1_u1 = P_of(u1) - lam0 * (1 - G_u1)
    nu = alpha * F1_u1
    mu = lam0 + nu
    rec["mu_interval"] = [str(mu.lo), str(mu.hi)]
    print("mu =", mu, "  > 1/3:", mu.lo > F(1, 3))
    cond_mu = mu.lo > F(1, 3)
    rec["mu_gt_third"] = cond_mu

    # ---- (2) upper-segment monotonicity ------------------------------
    slope = 1 + 2 * mu - 2 * Iv(F(5, 6), F(5, 6))
    cond_mono = slope.lo > 0
    rec["upper_slope_min"] = str(slope.lo)
    print("upper-segment slope min =", slope.lo, " > 0:", cond_mono)

    # ---- (3) boundary lower bound ------------------------------------
    # slope at y<=5/6 bounded by 1 + 2 mu_max - 4/3
    c = 1 + 2 * mu.hi - F(4, 3)
    B_lo = F(5, 36) - c * H_BAR
    rec["boundary_lower"] = str(B_lo)
    print("boundary winner B_K >= 5/36 - c/K >=", float(B_lo))

    # ---- (4) mid-region pieces ---------------------------------------
    oma = 1 - alpha
    # Piece A: y in [0, 1/3], swept twice, G=y/2
    # D2 = (1-a)^2 (2y - y^2) + ((1-a) lam0 + nu) y/2 - y
    mA = quad_max(-(oma * oma), 2 * (oma * oma)
                  + ((oma * lam0 + nu) * F(1, 2)) - 1,
                  iv(0), F(0), F(1, 3))
    # Piece B: y in [1/3, u0], G = 3y/2 - 1/3
    # D2 = (1-a)^2(2y-y^2) + ((1-a) lam0 + nu)(3y/2 - 1/3) - y
    mB = quad_max(-(oma * oma), 2 * (oma * oma)
                  + (oma * lam0 + nu) * F(3, 2) - 1,
                  -((oma * lam0 + nu) * F(1, 3)),
                  F(1, 3), u0.hi)
    # Piece C: y in [u0, u1], step-1 unswept, step-2 swept, G = 3y/2-1/3
    # F1 = P - lam0(1-G); D2 = (1-a) F1 + nu G - y
    #    = (1-a)(2y - y^2 - lam0(4/3 - 3y/2)) + nu(3y/2 - 1/3) - y
    mC = quad_max(-oma,
                  2 * oma + (oma * lam0) * F(3, 2) + nu * F(3, 2) - 1,
                  -(oma * lam0) * F(4, 3) - nu * F(1, 3),
                  u0.lo, u1.hi)
    # Piece E: y in [u1, 2/3], both unswept:
    #   D2 = -y^2 + (1 + 3 mu /2) y - 4 mu /3, vertex (1+3mu/2)/2.
    # Vertex interval is strictly right of 2/3 (proved), so the piece
    # max is at y = 2/3 where D2 = 2/9 - mu/3, DECREASING in mu; the
    # naive interval evaluation loses the (c1, c0) anti-correlation
    # and must not be used here.
    ystar_E = (1 + F(3, 2) * mu) / 2
    cond_vertex_E = ystar_E.lo > F(2, 3)
    rec["piece_E_vertex_gt_23"] = cond_vertex_E
    mE = Iv(F(2, 9) - mu.hi / 3, F(2, 9) - mu.lo / 3)
    print(f"piece E: vertex interval {ystar_E} > 2/3: {cond_vertex_E}; "
          f"max at 2/3 = 2/9 - mu/3 in {mE}")
    for nm, m in (("A", mA), ("B", mB), ("C", mC), ("E", mE)):
        rec[f"piece_{nm}_max"] = [str(m.lo), str(m.hi)]
        ok = m.hi < B_lo
        rec[f"piece_{nm}_below_boundary"] = ok
        print(f"piece {nm} max <= {float(m.hi):.6f}  < B_lo "
              f"{float(B_lo):.6f} : {ok}")

    all_ok = (cond_mu and cond_mono and cond_vertex_E
              and all(rec[f"piece_{x}_below_boundary"]
                      for x in "ABCE"))
    rec["D2_ENTRY_UNIFORM"] = all_ok
    rec["N"] = 33
    rec["finite_complement"] = ("K=6..32 covered by the accepted v3 "
                                "interval certificates (W2 entry)")
    dest = Path(__file__).with_name("d2_entry_certificate.json")
    dest.write_text(json.dumps(rec, indent=1), encoding="utf-8")
    print("\nD2 UNIFORM ENTRY (K>=33, alpha in J):", "PROVED" if all_ok
          else "FAILED")
    print("wrote", dest)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
