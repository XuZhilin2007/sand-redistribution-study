# -*- coding: utf-8 -*-
"""v04 all-K main line, step 1: EXACT continuum limit of the first two
updates and the W0 limit margins, as explicit polynomials in alpha.

Limit model (K -> infinity, canonical init, XO kernel, alpha fixed):
  density f0(y) = 2 - 2y
  step 1 selects y = 1/2 (F0* = 3/4)
  step 2 selects yM = 1/2 + 9*alpha/16 (prefix-lemma limit)
  respray at step t adds alpha*F_t * g_d(y) with the PIECEWISE CONSTANT
  density g_d = 1/2, 3/2, 2, 0 on [0,1/3), [1/3,2/3), [2/3,5/6), [5/6,1]
  (NOTE: density, not the CDF G)
  U_t = int_0^1 (f_t - 1)^2 dy

f_t is piecewise linear in y with alpha-polynomial coefficients and
breakpoints in {1/3, 1/2, yM(alpha), 2/3, 5/6} -- all integrals are
alpha-polynomials.  Positivity of the limit W0 margins on J is
certified with symkernel (exact).
"""
from __future__ import annotations

import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "v04_w_interval"))
from symkernel import pos_closed  # noqa: E402


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(n)]


def psub(a, b):
    return padd(a, [-c for c in b])


def pmul(a, b):
    if not a or not b:
        return []
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def peval(a, x):
    r = F(0)
    for c in reversed(a):
        r = r * x + c
    return r


def pscal(a, s):
    return [c * s for c in a]


def pl_int(pl):
    """int of value(y)dy over pieces (y0,y1 affine in alpha)."""
    tot = [F(0)]
    for y0, y1, sl, ic in pl:
        s2 = pscal(psub(pmul(y1, y1), pmul(y0, y0)), F(1, 2))
        tot = padd(tot, padd(pmul(sl, s2),
                             pmul(ic, psub(y1, y0))))
    return tot


def pl_sq_int(pl):
    """int of value(y)^2 dy over pieces."""
    tot = [F(0)]
    for y0, y1, sl, ic in pl:
        s2 = psub(pmul(y1, y1), pmul(y0, y0))
        s3 = psub(pmul(pmul(y1, y1), y1), pmul(pmul(y0, y0), y0))
        tot = padd(tot, padd(
            pscal(pmul(pmul(sl, sl), s3), F(1, 3)),
            padd(pmul(pmul(sl, ic), s2),
                 pmul(pmul(ic, ic), psub(y1, y0)))))
    return tot


A = [F(0), F(1)]
OMA = [F(1), F(-1)]
OMA2 = pmul(OMA, OMA)
HALF, YM, TWOTHIRD, FIVESIXTH = [F(1, 2)], padd([F(1, 2)],
                                                 pscal(A, F(9, 16))), \
    [F(2, 3)], [F(5, 6)]
F0STAR = F(3, 4)

# segment table: (lo, hi, g_d)
SEG = [(F(0), F(1, 3), F(1, 2)), (F(1, 3), F(1, 2), F(3, 2)),
       (F(1, 2), F(2, 3), F(3, 2)), (F(2, 3), F(5, 6), F(2)),
       (F(5, 6), F(1), F(0))]


def f1_pieces():
    """f1 = (1-a) f0 + a F0 g_d below 1/2 ; f0 + a F0 g_d above."""
    out = []
    for lo, hi, gd in SEG:
        m = OMA if hi <= F(1, 2) else [F(1)]
        out.append(([F(lo)], [F(hi)], pscal(m, F(-2)),
                    padd(pscal(m, F(2)), pscal(A, F0STAR * gd))))
    return out


def F1_poly(upper=False):
    f1 = f1_pieces()
    if not upper:
        # integrate [0,1/3] + [1/3,1/2] + [1/2, yM(1/2+9a/16)]
        ysel = padd([F(1, 2)], pscal(A, F(9, 16)))
        return padd(pl_int(f1[:2]),
                    pl_int([(HALF, ysel) + tuple(f1[2][2:])]))
    # upper: F1 = int over [0,1/3]+[1/3,1/2]+[1/2,2/3]+[2/3, yU]
    ysel = padd([F(1, 2)], pscal(A, F(3, 4)))
    return padd(pl_int(f1[:3]),
                pl_int([([F(2, 3)], ysel) + tuple(f1[3][2:])]))


def f2_pieces(upper=False):
    """f2: original density swept k(y) times; respray masses:
    step-1 respray swept again iff y <= ySel; step-2 respray fresh.
    upper=False: second selection yM = 1/2 + 9a/16 in the middle
    segment (valid for alpha <= 16/63 where the middle peak wins);
    upper=True: yU = 1/2 + 3a/4 in [2/3, 5/6] (the prefix-lemma y_U)."""
    F1 = F1_poly(upper)
    ysel = padd([F(1, 2)], pscal(A, F(3, 4))) if upper else \
        padd([F(1, 2)], pscal(A, F(9, 16)))
    out = []
    if not upper:
        for lo, hi, gd in SEG:
            if hi <= F(1, 2):
                m, mres = OMA2, OMA      # both steps sweep
                out.append(([F(lo)], [F(hi)], pscal(m, F(-2)),
                            padd(pscal(m, F(2)),
                                 padd(pmul(A, pscal(mres, F0STAR * gd)),
                                      pmul(pmul(A, F1), [gd])))))
            elif lo == F(1, 2):
                m, mres = OMA, OMA       # step-2 sweep; respray swept
                out.append((HALF, ysel, pscal(m, F(-2)),
                            padd(pscal(m, F(2)),
                                 padd(pmul(A, pscal(mres, F0STAR * gd)),
                                      pmul(pmul(A, F1), [gd])))))
                m, mres = [F(1)], [F(1)]  # unswept above
                out.append((ysel, TWOTHIRD, pscal(m, F(-2)),
                            padd(pscal(m, F(2)),
                                 padd(pmul(A, pscal(mres, F0STAR * gd)),
                                      pmul(pmul(A, F1), [gd])))))
            else:
                m, mres = [F(1)], [F(1)]
                out.append(([F(lo)], [F(hi)], pscal(m, F(-2)),
                            padd(pscal(m, F(2)),
                                 padd(pmul(A, pscal(mres, F0STAR * gd)),
                                      pmul(pmul(A, F1), [gd])))))
        return out
    # upper branch: ysel = 1/2 + 3a/4 in (2/3, 5/6)
    for lo, hi, gd in SEG:
        if hi <= F(1, 2):
            m, mres = OMA2, OMA
            out.append(([F(lo)], [F(hi)], pscal(m, F(-2)),
                        padd(pscal(m, F(2)),
                             padd(pmul(A, pscal(mres, F0STAR * gd)),
                                  pmul(pmul(A, F1), [gd])))))
        elif lo == F(1, 2) or lo == F(2, 3):
            hi_y = ysel if lo == F(2, 3) else TWOTHIRD
            m, mres = [F(1)], OMA if lo == F(1, 2) else OMA
            # (1/2,2/3]: original unswept by step 1, swept by step 2;
            # (2/3, ysel]: same
            out.append(([F(lo)], hi_y, pscal(OMA, F(-2)),
                        padd(pscal(OMA, F(2)),
                             padd(pmul(A, pscal(OMA, F0STAR * gd)),
                                  pmul(pmul(A, F1), [gd])))))
            if lo == F(2, 3):
                m2, mres2 = [F(1)], [F(1)]
                out.append((ysel, FIVESIXTH, pscal(m2, F(-2)),
                            padd(pscal(m2, F(2)),
                                 padd(pmul(A, pscal(mres2, F0STAR * gd)),
                                      pmul(pmul(A, F1), [gd])))))
        elif lo == F(5, 6):
            m, mres = [F(1)], [F(1)]
            out.append(([F(lo)], [F(hi)], pscal(m, F(-2)),
                        padd(pscal(m, F(2)),
                             padd(pmul(A, pscal(mres, F0STAR * gd)),
                                  pmul(pmul(A, F1), [gd])))))
    return out


def main():
    f2m = f2_pieces(upper=False)
    f2u = f2_pieces(upper=True)
    U0 = [F(1, 3)]
    f1 = f1_pieces()
    U1 = pl_sq_int([(y0, y1, sl, psub(ic, [F(1)]))
                    for y0, y1, sl, ic in f1])

    def U2_of(f2):
        return pl_sq_int([(y0, y1, sl, psub(ic, [F(1)]))
                          for y0, y1, sl, ic in f2])
    U2m, U2u = U2_of(f2m), U2_of(f2u)
    a0 = F(1, 4)
    print("middle-branch F1(1/4) =", float(peval(F1_poly(False), a0)),
          "(discrete 0.801025)")
    print("U2m(1/4) =", float(peval(U2m, a0)), "(discrete 0.150340)")
    # J and the LEFT candidate use the middle branch (alpha <= 16/63)
    lo, hi = F(1999, 8000), F(2001, 8000)
    print("[middle branch, valid a<=16/63]")
    print("U0-U2 > 0 on J:", pos_closed(psub(U0, U2m), lo, hi),
          " margins:", float(peval(psub(U0, U2m), lo)),
          float(peval(psub(U0, U2m), hi)))
    print("U1-U2 > 0 on J:", pos_closed(psub(U1, U2m), lo, hi),
          " margins:", float(peval(psub(U1, U2m), lo)),
          float(peval(psub(U1, U2m), hi)))
    for name, l2, h2 in [("left-candidate", F(59, 240), F(119, 480)),
                         ("up-to-16/63", F(1, 5), F(16, 63))]:
        print(f"  {name} [{float(l2):.6f},{float(h2):.6f}]: "
              f"U0-U2 {pos_closed(psub(U0, U2m), l2, h2)}, "
              f"U1-U2 {pos_closed(psub(U1, U2m), l2, h2)}")
    # the RIGHT candidate uses the UPPER branch (alpha >= 16/63)
    print("[upper branch, valid a>=16/63]")
    for name, l2, h2 in [("right-candidate", F(16, 63), F(7, 24)),
                         ("right-candidate-bracket", F(17, 60),
                          F(7, 24))]:
        print(f"  {name} [{float(l2):.6f},{float(h2):.6f}]: "
              f"U0-U2 {pos_closed(psub(U0, U2u), l2, h2)}, "
              f"U1-U2 {pos_closed(psub(U1, U2u), l2, h2)}")
    out = {"U0": ["1/3"], "U1": [str(c) for c in U1],
           "U2_middle": [str(c) for c in U2m],
           "U2_upper": [str(c) for c in U2u]}
    Path(__file__).with_name("continuum_U_polys.json").write_text(
        __import__("json").dumps(out, indent=1), encoding="utf-8")
    print("wrote continuum_U_polys.json")


if __name__ == "__main__":
    main()
