# -*- coding: utf-8 -*-
"""Research direction B (2026-10-19/20 handover; FIXED per owner review 7):
all-K, all-time stopping theorem for the LOCK class along the canonical
families K = 6n + r, r in {0,1,2}.

REVIEW-7 CORRECTIONS (this version):
* The old claim C5 (S'^2 + (U2-U_inf)E' <= 0) was certified in the WRONG
  direction (the generator proved >= 0) and is mathematically false --
  at K=12: S'^2 + (U2-U_inf)E' = 205209/1073741824 > 0, and the parabola
  vertex z* = -S'/E' = 226240/157449 > 1 lies OUTSIDE the reachable range
  z = (3/4)^{t-t_lock-1} in (0,1].  Replaced by C5': -(S'+E') > 0, which
  with E' > 0 gives Q'(z) = 2(S'+E'z) <= 2(S'+E') < 0 on (0,1], so the
  tail minimum is Q(1) and C4 (Q(1) >= U2-U_inf) alone completes the
  all-time bound.  Both owner anchor values are asserted at K=12.
* The second-selection derivation is now the two-peak lemma (D_1 has TWO
  interior peaks, at x = 41/64 and x = 11/16, heights 1425/4096 and
  1424/4096, difference exactly 1/4096 -- the band-3 piece is NOT
  decreasing).  The closed form (owner):
      s_1(K) = floor((41K+31)/64)  for even K (r = 0, 2),
      s_1(K) = floor((41K+26)/64)  for odd  K (r = 1),
  with M_0 = 3/16 constant for even K and M_0 = 3/16 - 1/(8K) - 1/(16K^2)
  for odd K.  r=2's M_0 is CONSTANT (s_0 = 3n+1 = K/2), so only r=1 has
  the n-dependent vertex; the rounding is absorbed exactly (see the lemma
  proof in docs/V03_B_LOCK_ALLK.md).
* Engine cell tables now handle the XO piecewise breakpoints that fall
  INSIDE bins: r=1 has mixed bins at x=1/3 (g = 7/(6K)) and x=2/3
  (g = 5/(3K)); r=2 has mixed bins at x=1/3 (g = 5/(6K)) and x=2/3
  (g = 11/(6K)).  r=0's clean band split does NOT carry over.
* Certificates are generated for all three residues: 96 sub-families
  (r in {0,1,2} x nu in 0..31, n = 32u + nu), each with claims
      C1: -S' > 0                 C4: Q(1) - (U2 - U_inf) >= 0
      C2: U(0)-U(2) >= 0          C5': -(S'+E') > 0
      C3: U(1)-U(2) >= 0          C6: E' > 0
      C3b (r=1,2): U(3)-U(2) >= 0
  Verification binds: expected sub-family set exactly (96, unique),
  claim-set exactness, finite-prefix completeness from the claimed start,
  and the generator source hash.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from fractions import Fraction as F
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "s4_w1_s1_xo_precision"))
sys.path.insert(0, str(REPO / "experiments" / "v03_p_plateau_allk"))

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(1000000)

from run_a import (RF, hybrid_certificate, pp_eval,  # noqa: E402
                   verify_certificate)
import run_a  # noqa: E402
import replay  # noqa: E402

run_a.CANCEL = True   # gcd-cancel + sign-normalize every RF

RESULTS = Path(__file__).resolve().parent / "results"
RHO = F(3, 4)
ALPHA = F(1, 4)
ZERO = RF([F(0)], [F(1)])
N_CHECK = 300


def rf(k):
    return RF([F(k)], [F(1)])


def rf_lin(a, b):
    return RF([F(b), F(a)], [F(1)])


def s_i(x):
    return x * (x + 1) / 2


def s_i2(x):
    return x * (x + 1) * (2 * x + 1) / 6


NBIG = 10 ** 6


def v(x):
    return pp_eval(x.num, NBIG) / pp_eval(x.den, NBIG)


class Cell:
    """Bins [lo, hi] (RF endpoints) with bin mass g (constant on the cell;
    mixed breakpoint bins are their own cells) and profile coefficients:
    p_i = A*p0_i + B*g_i (p0_i = 2/K-(2i-1)/K^2)."""

    __slots__ = ("lo", "hi", "g", "A", "B")

    def __init__(self, lo, hi, g, A, B):
        self.lo, self.hi, self.g, self.A, self.B = lo, hi, g, A, B

    def N(self):
        return self.hi - self.lo + 1

    def I(self):
        return s_i(self.hi) - s_i(self.lo) + self.lo

    def I2(self):
        return s_i2(self.hi) - s_i2(self.lo - 1)

    def sum_p0(self, K):
        return self.N() * 2 / K - (self.I() * 2 - self.N()) / (K * K)

    def sum_g(self):
        return self.g * self.N()

    def sum_p0sq(self, K):
        N, I, I2 = self.N(), self.I(), self.I2()
        return N * 4 / (K * K) - (8 * I - 4 * N) / (K ** 3) + \
            (4 * I2 - 4 * I + N) / (K ** 4)

    def sum_gsq(self):
        return self.g * self.g * self.N()

    def sum_p0g(self, K):
        return self.g * self.sum_p0(K)

    def sum_p(self, K):
        return self.A * self.sum_p0(K) + self.B * self.sum_g()

    def sum_psq(self, K):
        return (self.A * self.A * self.sum_p0sq(K) +
                self.A * self.B * 2 * self.sum_p0g(K) +
                self.B * self.B * self.sum_gsq())


def evolve(profile, s, M):
    """One sweep with selection s, respray mass M = alpha*F_s: prefix bins
    i <= s scale by 3/4; EVERY bin with g_i != 0 gains M*g_i."""
    out = []
    for c in profile:
        has_g = len(c.g.num) > 0
        Bn = c.B + M if has_g else c.B
        if v(c.hi) <= v(s):
            out.append(Cell(c.lo, c.hi, c.g, c.A * RHO,
                            c.B * RHO + M if has_g else c.B * RHO))
        elif v(c.lo) > v(s):
            out.append(Cell(c.lo, c.hi, c.g, c.A, Bn))
        else:
            out.append(Cell(c.lo, s, c.g, c.A * RHO,
                            c.B * RHO + M if has_g else c.B * RHO))
            out.append(Cell(s + 1, c.hi, c.g, c.A, Bn))
    return out


def sum_p_upto(profile, s, K):
    tot = ZERO
    for c in profile:
        if v(c.lo) > v(s):
            continue
        hi = s if v(c.hi) > v(s) else c.hi
        tot = tot + Cell(c.lo, hi, c.g, c.A, c.B).sum_p(K)
    return tot


def sum_psq_all(profile, K):
    tot = ZERO
    for c in profile:
        tot = tot + c.sum_psq(K)
    return tot


# ---------------- sub-family tables (r, nu) ----------------
def subfamilies():
    """All 96 sub-families: (r, nu, u0, sweeps, cells_spec).  n = 32u+nu;
    K = 192u + 6nu + r; m = 160u + 5nu + r; band edges 2n, 4n; sweeps
    use the owner closed form s_1 = floor((41K+c)/64), c = 31 (even K),
    26 (odd K).  Cells include the mixed breakpoint bins for r=1,2."""
    out = []
    for r in (0, 1, 2):
        for nu in range(32):
            n_min = {0: 2, 1: 4, 2: 1}[r]
            u0 = 0 if nu >= n_min else 1
            Kc = 192 * 0 + 6 * nu + r          # K at u=0 for offsets
            c31 = (41 * Kc + (31 if r in (0, 2) else 26)) // 64
            s0_b = 3 * nu + (0 if r in (0, 1) else 1)
            s1 = (123, c31)
            sweeps = [(96, s0_b), s1]
            if r in (1, 2):
                ell_b = 5 * nu + (0 if r == 1 else 1)
                sweeps.append((160, ell_b))
            out.append({"r": r, "nu": nu, "u0": u0, "sweeps": sweeps})
    return out


def build_cells(r, nu):
    """Base cells for sub-family (r, nu); variable u.  n = 32u+nu."""
    two_n = rf_lin(64, 2 * nu)          # 2n
    four_n = rf_lin(128, 4 * nu)        # 4n
    K = rf_lin(192, 6 * nu + r)
    m = rf_lin(160, 5 * nu + r)
    q = rf(6 - r) / (K * 3)
    A0, B0 = rf(1), rf(0)
    g1 = F(1, 2) / K
    g2 = F(3, 2) / K
    g3 = F(2) / K
    if r == 0:
        cells = [
            Cell(rf(1), two_n, g1, A0, B0),
            Cell(two_n + 1, four_n, g2, A0, B0),
            Cell(four_n + 1, m - 1, g3, A0, B0),
        ]
    elif r == 1:
        # mixed bins: x=1/3 inside bin 2n+1 (g=7/(6K)); x=2/3 inside bin
        # 4n+1 (g=5/(3K))
        cells = [
            Cell(rf(1), two_n, g1, A0, B0),
            Cell(two_n + 1, two_n + 1, F(7) / (6 * K), A0, B0),
            Cell(two_n + 2, four_n, g2, A0, B0),
            Cell(four_n + 1, four_n + 1, F(5) / (3 * K), A0, B0),
            Cell(four_n + 2, m - 1, g3, A0, B0),
        ]
    else:
        # mixed bins: x=1/3 inside bin 2n+1 (g=5/(6K)); x=2/3 inside bin
        # 4n+2 (g=11/(6K)); note the pure band-2 cell ends at 4n+1
        cells = [
            Cell(rf(1), two_n, g1, A0, B0),
            Cell(two_n + 1, two_n + 1, F(5) / (6 * K), A0, B0),
            Cell(two_n + 2, four_n + 1, g2, A0, B0),
            Cell(four_n + 2, four_n + 2, F(11) / (6 * K), A0, B0),
            Cell(four_n + 3, m - 1, g3, A0, B0),
        ]
    cells.append(Cell(m, m, q, A0, B0))
    cells.append(Cell(m + 1, K, ZERO, A0, B0))
    return cells, K, m


def build_family(sub):
    r, nu, u0 = sub["r"], sub["nu"], sub["u0"]
    t_lock = 2 if r == 0 else 3
    cells, K, m = build_cells(r, nu)
    b = m / K
    c = b * 2 - b * b
    q = rf(6 - r) / (K * 3)
    L = c * q
    profiles = [cells]
    for (a, bb) in sub["sweeps"]:
        s = rf_lin(a, bb)
        Fs = sum_p_upto(profiles[-1], s, K)
        profiles.append(evolve(profiles[-1], s, Fs * ALPHA))
    pA = evolve(profiles[-1], m, c * ALPHA)
    pA_int = [cl for cl in pA if v(cl.hi) < v(m)]
    e_int = [Cell(cl.lo, cl.hi, cl.g, cl.A, cl.B - c) for cl in pA_int]
    mc = next(cl for cl in profiles[-1] if v(cl.lo) == v(m))
    x = mc.A * (2 / K - (m * 2 - 1) / (K * K)) + mc.B * q
    e_m = RHO * (x - L)
    # boundary-coordinate crossings for the real-orbit prefix (review 8):
    # x_t = p_{t,m} from the m-cell of each profile (bin m is resprayed at
    # every step, scaled only by the m/A step)
    p0m = 2 / K - (m * 2 - 1) / (K * K)

    def mcell_x(profile):
        cl = next(cc for cc in profile if v(cc.lo) == v(m))
        return cl.A * p0m + cl.B * q

    x2 = mcell_x(profiles[2])
    x3 = mcell_x(profiles[3]) if t_lock == 3 else None
    U0 = K * sum_psq_all(profiles[0], K) - 1
    U1 = K * sum_psq_all(profiles[1], K) - 1
    U2 = K * sum_psq_all(profiles[2], K) - 1
    U3 = K * sum_psq_all(profiles[3], K) - 1 if t_lock == 3 else None
    g2s, suf = ZERO, ZERO
    for cl in profiles[0]:
        if v(cl.hi) < v(m):
            g2s = g2s + cl.sum_gsq()
        if v(cl.lo) > v(m):
            suf = suf + cl.sum_p0sq(K)
    Uinf = K * (c * c * g2s + L * L + suf) - 1
    Ssum, Esum = ZERO, ZERO
    for cl in e_int:
        Ssum = Ssum + cl.A * c * cl.sum_p0g(K) + cl.B * c * cl.sum_gsq() \
            - (cl.A * cl.sum_p0(K) + cl.B * cl.sum_g()) / K
        Esum = Esum + cl.sum_psq(K)
    Sp = K * (Ssum + e_m * (L - 1 / K))
    Ep = K * (Esum + e_m * e_m)
    return {"K": K, "m": m, "c": c, "q": q, "L": L, "u0": u0,
            "t_lock": t_lock, "U0": U0, "U1": U1, "U2": U2, "U3": U3,
            "Uinf": Uinf, "Sp": Sp, "Ep": Ep, "x2": x2, "x3": x3}


def claims_for(D):
    U2, Sp, Ep, Uinf = D["U2"], D["Sp"], D["Ep"], D["Uinf"]
    K, h, L = D["K"], 1 / D["K"], D["L"]
    claims = {
        "C1_-Sp>0": (Sp * (-1), True),
        "C2_U0-U2>=0": (D["U0"] - U2, False),
        "C3_U1-U2>=0": (D["U1"] - U2, False),
        "C4_Q1>=U2-Uinf": (Sp * 2 + Ep - (U2 - Uinf), False),
        "C5p_-(Sp+Ep)>0": ((Sp + Ep) * (-1), True),
        "C6_Ep>0": (Ep, True),
        # real-orbit prefix crossings (review 8): selection at t is m iff
        # x_t > h (T1/T4 scalar closure); persistence thereafter is T3
        "P1_L-h>0": (L - h, True),
    }
    if D["t_lock"] == 2:
        claims["P2_x2-h>0"] = (D["x2"] - h, True)       # r=0: s_2 = m
    else:
        claims["P2_h-x2>0"] = (h - D["x2"], True)       # r=1,2: s_2 = ell
        claims["P3_x3-h>0"] = (D["x3"] - h, True)       # then s_3 = m
    if D["U3"] is not None:
        claims["C3b_U3-U2>=0"] = (D["U3"] - U2, False)
    return claims


# ---------------- canonical prefix ----------------
def s1_closed(K):
    """Owner closed form: floor((41K+31)/64) even K, floor((41K+26)/64)
    odd K."""
    return (41 * K + 31) // 64 if K % 2 == 0 else (41 * K + 26) // 64


def expected_prefix(K):
    """canonical selection prefix (5 steps) for LOCK K >= 7."""
    n, r = divmod(K, 6)
    m = (5 * K + 5) // 6
    s0 = 3 * n + (1 if r == 2 else 0)
    s1 = s1_closed(K)
    if r == 1 and n <= 3:
        return [s0, s1, m, m, m]          # K in {7,13,19}: no ell step
    if r == 0:
        return [s0, s1, m, m, m]
    return [s0, s1, 5 * n + (0 if r == 1 else 1), m, m]


def lemma_t2():
    """All-K boundary-entry lemma (owner review-8 analytic patch; numeric
    statements corrected per review 9), EXACT rational arithmetic, valid
    for every LOCK K >= 33: the t=2 grid argmax of D_2 lies in {l, m},
    hence t* = 2 for the REAL orbit (K in [7, 33) is covered by the
    committed replay check).  Chain:
      (a) D_1 <= 1/6 everywhere.  Peak heights are NOT uniform in K:
          even K (M_0 = 3/16) peaks are 657/4096 and 656/4096; odd K
          has M_0 in [M(37), 3/16) and f1, f2 DECREASE in M, so the
          uniform bound is max(f1(M37), f2(M37), f1(3/16)) with
          f1(M37) = 1212073/7496644 < 1/6 (owner-checked; K=37 grid max
          8184/50653 > 657/4096 anchors that 657/4096 is not universal).
          Non-peak pieces are uniformly <= 9/64: swept band 2 attains
          9/64 at (M, x) = (3/16, 1/2); swept band 1 <= 11/96 (attained
          at K=36, x=1/3) < 9/64; suffix and bin m <= 5/36 < 9/64.
      (b) j <= s_1: D_2,j <= (3/4)D_1,j + (j/K)(M_1 - 1/4) <= 1/8
          (G <= id on x <= 2/3, s_1/K < 2/3, M_1 <= 1/4);
      (c) anchor D_2,m = D_0,m = b(1-b) > 1/8  (b in [5/6, 16/19];
          worst r=2, K=38);
      (d) j > s_1: D_2 = -x^2+(1+3S/2)x-4S/3 (band 2), -x^2+(1+2S)x-5S/3
          (band 3) with S = M_0+M_1 > 1/3: band-2 vertex > 1/2+3/12 =
          3/4 > 2/3 and band-3 vertex >= 1/2+1/3 = 5/6, so both are
          increasing on (s_1/K, 5/6]; suffix = x(1-x) decreasing ->
          grid max over j > s_1 in {l, m}.
    Returns failure list (empty = verified)."""
    fails = []
    sixth = F(1, 6)
    # M_0 range for K >= 33: even K M_0 = 3/16; odd K >= 37 increasing in
    # K, so M_0 in [M37, 3/16]
    M37 = F(3, 16) - F(1, 8 * 37) - F(1, 16 * 37 * 37)
    # (a) true peak heights and their monotonicity (f1, f2 decrease in M,
    # so the uniform peak bound on the M-range is attained at M = M37 for
    # odd K and at M = 3/16 for even K)
    f1 = lambda M: (F(1, 2) + 3 * M / 4) ** 2 - 4 * M / 3
    f2 = lambda M: (F(1, 2) + M) ** 2 - 5 * M / 3
    d1 = lambda M: F(9, 8) * M - F(7, 12)      # f1'
    d2 = lambda M: 2 * M - F(2, 3)             # f2'
    if not (d1(M37) < 0 and d1(F(3, 16)) < 0):
        fails.append("f1 not decreasing on M-range")
    if not (d2(M37) < 0 and d2(F(3, 16)) < 0):
        fails.append("f2 not decreasing on M-range")
    if not (f1(M37) == F(1212073, 7496644) and f1(M37) < sixth
            and f1(F(3, 16)) == F(657, 4096) and f1(F(3, 16)) < sixth):
        fails.append("f1 peak bound wrong or >= 1/6")
    if not (f2(M37) < sixth and f2(F(3, 16)) == F(656, 4096)
            and f2(F(3, 16)) < sixth):
        fails.append("f2 peak bound wrong or >= 1/6")
    # peak gap (Lemma S2 dominance): Delta decreasing on [33/196, 3/16],
    # minimum at the right endpoint = exactly 1/4096
    dd = lambda M: F(1, 12) * M - F(7, 16) * M * M
    ddd = lambda M: F(1, 12) - F(7, 8) * M
    if not (ddd(F(33, 196)) < 0 and ddd(F(3, 16)) < 0):
        fails.append("Delta not decreasing on [33/196, 3/16]")
    if not (dd(F(33, 196)) >= F(1, 4096) and dd(F(3, 16)) == F(1, 4096)):
        fails.append("peak gap < 1/4096")
    # non-peak pieces: uniform bound D_1 <= 9/64 = 576/4096
    #   swept band 2: 1/16 + 5M/12, increasing in M, max 9/64 at 3/16
    #   swept band 1: -(3/4)x^2 + (1/2 + M/2)x <= 11/96 at (3/16, 1/3)
    #     (K=36 attains it), and 11/96 < 9/64
    #   suffix and bin m: x(1-x) <= 5/36 < 9/64 (b >= 5/6 there)
    if not (F(1, 16) + 5 * F(3, 16) / 12 == F(9, 64)):
        fails.append("swept band-2 max != 9/64")
    if not (-F(3, 4) * F(1, 9) + (F(1, 2) + F(3, 32)) / 3 == F(11, 96)
            and F(11, 96) < F(9, 64)):
        fails.append("swept band-1 max != 11/96 or >= 9/64")
    if not (F(5, 36) < F(9, 64) and F(9, 64) < sixth):
        fails.append("non-peak uniform bound 9/64 fails")
    # (b) s_1/K < 2/3 for K >= 33 (so G <= id applies below s_1)
    if not (F(41, 64) + F(31, 64 * 33) < F(2, 3)):
        fails.append("s_1/K >= 2/3")
    # (c) anchor b(1-b) > 1/8 for K >= 33: b = m/K with
    #     r=0: b = 5/6; r=1: b = (5n+1)/(6n+1) decreasing in n, max at
    #     K=37 -> 31/37; r=2: (5n+2)/(6n+2) decreasing, max at K=38 ->
    #     16/19; overall b in [5/6, 16/19], and b(1-b) decreasing in b
    blo, bhi = F(5, 6), F(16, 19)
    if not (F(31, 37) < bhi and min(blo * (1 - blo), bhi * (1 - bhi))
            > F(1, 8)):
        fails.append("anchor D_2,m <= 1/8")
    # (d) S = M_0 + M_1 > 1/3 with M_1 >= (F_{0,s1} - F_{0,s0}/4)/4,
    #     F_{0,s1} >= near_cdf(5/8) = 55/64 (s1/K >= 5/8 for K >= 33),
    #     F_{0,s0} <= 3/4, M_0 >= M37
    M1lb = (F(55, 64) - F(3, 16)) / 4
    if not (F(5, 8) <= F(41, 64) - F(33, 64 * 33)):
        fails.append("s1/K lower bound < 5/8")
    if not (M37 + M1lb > F(1, 3)):
        fails.append("M_0 + M_1 <= 1/3")
    # with S > 1/3: band-2 vertex > 1/2 + 3/4*(1/3) = 3/4 > 2/3 (review-9
    # correction: 3/4, not 7/8); band-3 vertex >= 1/2 + 1/3 = 5/6 (right
    # end of band 3) -> both increasing
    if not (F(1, 2) + 3 * F(1, 3) / 4 == F(3, 4) and F(3, 4) > F(2, 3)
            and F(1, 2) + F(1, 3) >= F(5, 6)):
        fails.append("vertex bounds fail")
    return fails


def verify_prefix_pattern(nmax=N_CHECK):
    """Replay-verify the canonical selection prefix for every LOCK K in
    [7, 6*nmax]."""
    for K in range(7, 6 * nmax + 1):
        if K % 6 not in (0, 1, 2):
            continue
        P = replay.exact_bin_probs(replay.near_cdf_exact, K)
        g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
        p = P
        expect = expected_prefix(K)
        for t in range(5):
            D = replay.cumsum_excess(p, K)
            s = replay.argmax_first(D)
            if s != expect[t]:
                return (K, t, s, expect[t])
            p = replay.update_state(p, g, s, ALPHA)
    return None


def numeric_reference(K):
    """Exact per-K values for anchor comparison."""
    P = replay.exact_bin_probs(replay.near_cdf_exact, K)
    g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
    m = replay.anchor_m_b(K)[0]
    states = [P]
    sels = []
    for t in range(6):
        D = replay.cumsum_excess(states[-1], K)
        s = replay.argmax_first(D)
        sels.append(s)
        states.append(replay.update_state(states[-1], g, s, ALPHA))
    U = [K * sum(y * y for y in st) - 1 for st in states[:4]]
    x = [st[m - 1] for st in states[:4]]
    t_lock = next(t for t in range(5)
                  if all(sels[u] == m for u in range(t, 5)))
    b = F(m, K)
    c = 2 * b - b * b
    pinf = [c * gi for gi in g[:m - 1]] + [c * g[m - 1]] + list(P[m:])
    e = [states[t_lock + 1][i] - pinf[i] for i in range(K)]
    return {"sels": sels, "m": m, "U": U, "x": x, "t_lock": t_lock,
            "Sp": K * sum(e[i] * (pinf[i] - F(1, K)) for i in range(K)),
            "Ep": K * sum(y * y for y in e),
            "Uinf": K * sum(y * y for y in pinf) - 1}


# ---------------- certificates ----------------
def certify(rf_, n0, strict):
    """Certify rf_ > 0 (strict) or >= 0 on n >= n0."""
    den = rf_.den
    if hybrid_certificate(den, n0, True):
        target, ds = rf_.num, 1
    elif hybrid_certificate(den * (-1), n0, True):
        target, ds = [-x for x in rf_.num], -1
    else:
        return None
    cert_num = hybrid_certificate(target, n0, strict)
    if cert_num is None:
        return None
    return {"num": cert_num, "den_sign": ds,
            "den_cert": hybrid_certificate(den, n0, True) if ds == 1
            else hybrid_certificate(den * (-1), n0, True)}


def verify_claim(rf_, stored, n0, strict):
    """Re-verify with the THEOREM's start and direction supplied by the
    caller (owner review 8): the stored certificate must agree with both,
    so a self-consistent certificate whose start was shifted later (with
    recomputed coefficients and emptied prefix) is rejected."""
    fails = []
    den = rf_.den if stored["den_sign"] == 1 else rf_.den * (-1)
    fails += [f"den: {x}" for x in verify_certificate(
        den, stored["den_cert"], expected_start=n0, expected_strict=True)]
    num = rf_.num if stored["den_sign"] == 1 else [-x for x in rf_.num]
    fails += [f"num: {x}" for x in verify_certificate(
        num, stored["num"], expected_start=n0, expected_strict=strict)]
    return fails


def _source_sha256():
    import hashlib
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def main():
    t0 = time.time()
    print("== B: LOCK all-K all-time stopping certificates (r=0,1,2; "
          "review-7 fixed) ==")
    # owner anchor: the OLD C5 was certified in the wrong direction and is
    # false at K=12 -- assert the exact review values
    D12 = build_family(next(s for s in subfamilies()
                            if s["r"] == 0 and s["nu"] == 2))
    u = 0
    Sp12 = pp_eval(D12["Sp"].num, u) / pp_eval(D12["Sp"].den, u)
    Ep12 = pp_eval(D12["Ep"].num, u) / pp_eval(D12["Ep"].den, u)
    U212 = pp_eval(D12["U2"].num, u) / pp_eval(D12["U2"].den, u)
    Ui12 = pp_eval(D12["Uinf"].num, u) / pp_eval(D12["Uinf"].den, u)
    assert Sp12 * Sp12 + (U212 - Ui12) * Ep12 == F(205209, 1073741824), \
        "owner anchor: old C5 value at K=12"
    assert -Sp12 / Ep12 == F(226240, 157449) and -Sp12 / Ep12 > 1, \
        "owner anchor: vertex outside (0,1] at K=12"
    assert -(Sp12 + Ep12) > 0, "owner anchor: repair claim at K=12"
    print("[B] owner anchors at K=12 reproduced exactly (old C5 false, "
          "vertex > 1, -(S'+E') > 0)")
    bad = verify_prefix_pattern()
    if bad is not None:
        raise SystemExit(f"[B] ABORT: prefix pattern violated at K={bad[0]}"
                         f" step {bad[1]}: got {bad[2]}, expected {bad[3]}")
    print(f"[B] canonical prefix (s0, s1-closed-form, [ell,] m, ...) "
          f"replay-verified for every LOCK K in [7, {6 * N_CHECK}]")
    lt2 = lemma_t2()
    if lt2:
        raise SystemExit(f"[B] ABORT: boundary-entry lemma failed: {lt2}")
    print("[B] lemma t*=2 verified (K>=33 analytic: D_1<=1/6, swept<=1/8 "
          "< anchor b(1-b), M_0+M_1>1/3 monotone -> argmax in {l,m}; "
          "K<33 by replay)")
    out = {"families": [], "config": {
        "experiment": "v03_b_lock_allk",
        "prefix_check_K_max": 6 * N_CHECK}}
    subs = subfamilies()
    for sub in subs:
        r, nu, u0 = sub["r"], sub["nu"], sub["u0"]
        D = build_family(sub)
        for u_ in (u0, u0 + 1, u0 + 7):
            n = 32 * u_ + nu
            K = 6 * n + r
            ref = numeric_reference(K)
            if ref["t_lock"] != D["t_lock"]:
                raise SystemExit(f"[B] ABORT r={r} nu={nu} K={K}: t_lock "
                                 f"{ref['t_lock']} != {D['t_lock']}")
            checks = [("U0", "U", 0), ("U1", "U", 1), ("U2", "U", 2),
                      ("Sp", "Sp", None), ("Ep", "Ep", None),
                      ("Uinf", "Uinf", None), ("x2", "x", 2)]
            if D["U3"] is not None:
                checks += [("U3", "U", 3), ("x3", "x", 3)]
            for key, grp, idx in checks:
                symv = D[key]
                val = pp_eval(symv.num, u_) / pp_eval(symv.den, u_)
                numv = ref[grp][idx] if idx is not None else ref[grp]
                if val != numv:
                    raise SystemExit(f"[B] ABORT r={r} nu={nu} K={K}: "
                                     f"{key} symbolic {val} != {numv}")
        certs = {}
        for name, (rf_, strict) in claims_for(D).items():
            st = certify(rf_, u0, strict)
            if st is None:
                raise SystemExit(f"[B] ABORT r={r} nu={nu} {name}: no "
                                 f"certificate (num degree "
                                 f"{len(rf_.num) - 1})")
            vf = verify_claim(rf_, st, u0, strict)
            if vf:
                raise SystemExit(f"[B] ABORT r={r} nu={nu} {name}: {vf}")
            certs[name] = {"num": st["num"], "den_sign": st["den_sign"],
                           "den_cert": st["den_cert"]}
        out["families"].append({
            "r": r, "nu": nu, "u0": u0, "t_lock": D["t_lock"],
            "claims": certs,
            "num_degrees": {k: len(p.num) - 1
                            for k, (p, _) in claims_for(D).items()},
        })
        if len(out["families"]) % 16 == 0:
            print(f"[B] {len(out['families'])}/96 sub-families certified")
    out["git"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                                capture_output=True, text=True).stdout.strip()
    out["source_sha256"] = _source_sha256()
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "b_lock_cert.json").write_text(json.dumps(out, indent=1),
                                              encoding="utf-8")
    fails = verify_from_json()
    if fails:
        raise SystemExit(f"[B] ABORT: JSON re-verification failed: "
                         f"{fails[:5]}")
    print(f"[B] PASS: 96/96 sub-families certified + JSON re-verified "
          f"({time.time() - t0:.1f}s)")
    return 0


def verify_from_json():
    """Hardened re-verification (owner review 7): exact expected
    sub-family set (r, nu) for r in {0,1,2} x nu 0..31 with the correct
    u0, unique; claim sets exact per t_lock; every stored certificate
    re-verified with full finite-prefix coverage; source hash bound."""
    data = json.loads((RESULTS / "b_lock_cert.json").read_text(encoding="utf-8"))
    fails = []
    if data.get("source_sha256") != _source_sha256():
        fails.append("source_sha256 mismatch: certificate file generated "
                     "by a different version of run_b2.py")
    expected = {(s["r"], s["nu"], s["u0"]) for s in subfamilies()}
    got = [(f["r"], f["nu"], f["u0"]) for f in data["families"]]
    if sorted(got) != sorted(expected):
        misses = sorted(expected - set(got))
        extras = sorted(set(got) - expected)
        fails.append(f"sub-family set mismatch: missing {misses[:6]}, "
                     f"extra/duplicated {extras[:6]}")
    for entry in data["families"]:
        key = (entry["r"], entry["nu"])
        sub = next((s for s in subfamilies()
                    if (s["r"], s["nu"]) == key), None)
        if sub is None or entry["u0"] != sub["u0"]:
            continue                      # already reported
        D = build_family(sub)
        claims = claims_for(D)
        stored = set(entry["claims"])
        if stored != set(claims):
            fails.append(f"r={entry['r']} nu={entry['nu']}: claim set "
                         f"{stored} != {set(claims)}")
            continue
        for name, (rf_, strict) in claims.items():
            st = entry["claims"][name]
            stored_cert = {
                "strict": st["num"]["strict"], "N": st["num"]["N"],
                "covers_from": st["num"].get("covers_from"),
                "shift_coeffs": st["num"]["shift_coeffs"],
                "finite": st["num"]["finite"]}
            den = rf_.den if st["den_sign"] == 1 else rf_.den * (-1)
            stored_den = {
                "strict": True, "N": st["den_cert"]["N"],
                "covers_from": st["den_cert"].get("covers_from"),
                "shift_coeffs": st["den_cert"]["shift_coeffs"],
                "finite": st["den_cert"]["finite"]}
            fails += [f"r={entry['r']} nu={entry['nu']} {name}: {m}"
                      for m in verify_certificate(
                          den, stored_den, expected_start=sub["u0"],
                          expected_strict=True)]
            num = rf_.num if st["den_sign"] == 1 else \
                [-x for x in rf_.num]
            fails += [f"r={entry['r']} nu={entry['nu']} {name}: {m}"
                      for m in verify_certificate(
                          num, stored_cert, expected_start=sub["u0"],
                          expected_strict=strict)]
    return fails


if __name__ == "__main__":
    raise SystemExit(main())
