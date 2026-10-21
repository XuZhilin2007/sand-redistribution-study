# -*- coding: utf-8 -*-
"""CYCLE-side all-K all-time optimal stopping (2026-10-20 sprint; review-10
REPAIRED: cross-term factor 2 in the quadratic identity, square-completion
certificate route, finite exceptions via the committed D certificates).

THEOREM E (target).  For canonical XO, alpha = 1/4, exact arithmetic,
min-index tie rule: every CYCLE class K satisfies
    U(p_t) >= U(p_2)   for all integers t >= 0.

REVIEW-10 CORRECTION (owner-caught core algebra error): the strengthened
decomposition is correct, but the quadratic identity carries a factor 2
on the cross term:
    U(p_t) = Q(x_t) + 2 z Phi(x_t) + K E2 z^2,   z = (3/4)^{t-2},
    Phi(x) = K gamma(x) GE - E1.
The previous version wrote Q + z Phi + K E2 z^2 and derived a false
lower bound R = Q + (3/4)Phi - U2 (owner counterexample K = 9, t = 3:
(U3-U2) - R(x3) = -12642529/241864704 < 0).  Everything downstream of
that step is replaced by the SQUARE-COMPLETION route (owner's patch,
pre-checked):
    U - U2 = A x^2 + B x + (C - U2) + 2 z (d x + f) + E z^2
           = A (x + (B + 2 d z)/(2 A))^2 + H(z),
    H(z) = (E - d^2/A) z^2 + (2 f - B d / A) z + (C - U2 - B^2/(4A)),
with A = [x^2]Q, B = [x]Q, C = Q(0), d = [x]Phi, f = Phi(0), E = K E2.
Claims S1-S4 (all RFs in u, shift certificates bound to caller
start/strictness):
    S1: A > 0;   S2: H2 := E - d^2/A > 0;
    S3: -H'(3/4) > 0   (H' = 2 H2 z + H1, H1 = 2 f - B d / A);
    S4: H(3/4) >= 0.
S1 makes the square >= 0 for EVERY real x (no reachable interval
needed); S2 + S3 make H strictly decreasing on [0, 3/4], so for
t >= 3 (z <= 3/4): U - U2 >= H(z) >= H(3/4) >= 0.  t = 2: z = 1 and
U(p_2) = Q(x_2) + 2 Phi(x_2) + K E2 EXACTLY (identity, anchor-checked).
t in {0, 1}: claims C0a/C0b.  Direct back-substitution anchors
U(p_t) = Q(x_t) + 2 z Phi(x_t) + K E2 z^2 are asserted at t = 2 and
t = 3 for every sub-family (both branches occur across residues).

Proof chain L1-L4 unchanged (each link exact):
  L1 prefix: s_0 = floor(K/2); s_1 = floor((41K+31)/64) even /
      floor((41K+26)/64) odd for K >= 33 (two-peak argument is
      parity-based, residue-free); K < 33 by replay.
  L2 entry: lemma_t2_cycle -- t* = 2 for CYCLE K >= 33 with anchor
      D_{2,m-1}; the odd-K M_0 lower bound is M(33) = 200/1089 (K = 33,
      35 are CYCLE; the LOCK lemma's M(37) bound does NOT carry over --
      review-10 correction), f1(M33) = 255907/1581228 < 1/6.
  L3 residence: registered theorem B part 3 / T3 canonical case;
      replay-verified.
  L4 decomp (STRENGTHENED D-decomp; two-branch residence only, no
      locked word): gamma(A(x)) = (3/4)gamma(x) + c/4 and
      gamma(B(x)) = (1-q/4)gamma(x); the B-branch respray (c-x)/4
      exactly cancels the gamma jump:
          p_{t,i} = gamma(x_t) g_i + (3/4)^{t-2} e_i   (i < m, t >= 2).

Finite difference set (below the certified family starts): K in
{9, 10, 11, 16, 17, 23}.  These are covered by the COMMITTED D
certificates (all-time: finite prefix + word license + infinite tail);
run_e re-verifies exactly those records through the INDEPENDENT D
verifier (verify_d_independent.py, second implementation) restricted to
the subset.  The previous 600-step window check is REMOVED (it was not
an all-time argument -- review-10 correction).
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

run_a.CANCEL = True

RESULTS = Path(__file__).resolve().parent / "results"
RHO = F(3, 4)
ALPHA = F(1, 4)
ZERO = RF([F(0)], [F(1)])
N_CHECK = 300
FINITE_EXCEPTIONS = [9, 10, 11, 16, 17, 23]


def rf(k):
    return RF([F(k)], [F(1)])


def rf_frac(x):
    return RF([F(x)], [F(1)])


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


# ---------------- sub-family tables ----------------
def s1_closed(K):
    return (41 * K + 31) // 64 if K % 2 == 0 else (41 * K + 26) // 64


def subfamilies():
    """(r, nu, u0, sweeps) with n = 32u+nu; square-completion route starts:
    r=3 n>=2 (K>=15), r=4 n>=3 (K>=22), r=5 n>=4 (K>=29); below that the
    finite difference set {9,10,11,16,17,23} is covered by the D
    certificates."""
    out = []
    for r in (3, 4, 5):
        n_min = {3: 2, 4: 3, 5: 4}[r]
        for nu in range(32):
            u0 = 0 if nu >= n_min else 1
            Kc = 6 * nu + r
            c31 = (41 * Kc + (31 if Kc % 2 == 0 else 26)) // 64
            s0_b = 3 * nu + (1 if r == 3 else 2)
            out.append({"r": r, "nu": nu, "u0": u0,
                        "sweeps": [(96, s0_b), (123, c31)]})
    return out


def build_cells(r, nu):
    two_n = rf_lin(64, 2 * nu)
    four_n = rf_lin(128, 4 * nu)
    K = rf_lin(192, 6 * nu + r)
    m = rf_lin(160, 5 * nu + r)
    q = rf(6 - r) / (K * 3)
    A0, B0 = rf(1), rf(0)
    g1 = F(1, 2) / K
    g2 = F(3, 2) / K
    g3 = F(2) / K
    if r == 3:
        cells = [
            Cell(rf(1), two_n + 1, g1, A0, B0),
            Cell(two_n + 2, four_n + 2, g2, A0, B0),
            Cell(four_n + 3, m - 1, g3, A0, B0),
        ]
    elif r == 4:
        cells = [
            Cell(rf(1), two_n + 1, g1, A0, B0),
            Cell(two_n + 2, two_n + 2, F(7) / (6 * K), A0, B0),
            Cell(two_n + 3, four_n + 2, g2, A0, B0),
            Cell(four_n + 3, four_n + 3, F(5) / (3 * K), A0, B0),
            Cell(four_n + 4, m - 1, g3, A0, B0),
        ]
    else:
        cells = [
            Cell(rf(1), two_n + 1, g1, A0, B0),
            Cell(two_n + 2, two_n + 2, F(5) / (6 * K), A0, B0),
            Cell(two_n + 3, four_n + 3, g2, A0, B0),
            Cell(four_n + 4, four_n + 4, F(11) / (6 * K), A0, B0),
            Cell(four_n + 5, m - 1, g3, A0, B0),
        ]
    cells.append(Cell(m, m, q, A0, B0))
    cells.append(Cell(m + 1, K, ZERO, A0, B0))
    return cells, K, m


def build_family(sub):
    """RFs in u for one sub-family: prefix U's, Q/Phi coefficient RFs,
    E2, x2 (for identity anchors)."""
    r, nu, u0 = sub["r"], sub["nu"], sub["u0"]
    cells, K, m = build_cells(r, nu)
    b = m / K
    c = b * 2 - b * b
    q = rf(6 - r) / (K * 3)
    profiles = [cells]
    for (a, bb) in sub["sweeps"]:
        s = rf_lin(a, bb)
        Fs = sum_p_upto(profiles[-1], s, K)
        profiles.append(evolve(profiles[-1], s, Fs * ALPHA))
    p2_int = [cl for cl in profiles[2] if v(cl.hi) < v(m)]
    p0m = 2 / K - (m * 2 - 1) / (K * K)

    def mcell_x(profile):
        cl = next(cc for cc in profile if v(cc.lo) == v(m))
        return cl.A * p0m + cl.B * q

    x2 = mcell_x(profiles[2])
    U0 = K * sum_psq_all(profiles[0], K) - 1
    U1 = K * sum_psq_all(profiles[1], K) - 1
    U2 = K * sum_psq_all(profiles[2], K) - 1
    gam2 = (c - x2) / (rf(1) - q)
    G1 = G2 = GE = E1 = E2 = ZERO
    for cl in p2_int:
        sg = cl.sum_g()
        sg2 = cl.sum_gsq()
        sp0 = cl.sum_p0(K)
        sp0g = cl.sum_p0g(K)
        sp02 = cl.sum_p0sq(K)
        Ae, Be = cl.A, cl.B - gam2
        G1 = G1 + sg
        G2 = G2 + sg2
        E1 = E1 + Ae * sp0 + Be * sg
        E2 = E2 + Ae * Ae * sp02 + 2 * Ae * Be * sp0g + Be * Be * sg2
        GE = GE + Ae * sp0g + Be * sg2
    Sfx = ZERO
    for cl in profiles[0]:
        if v(cl.lo) > v(m):
            Sfx = Sfx + (cl.sum_p0sq(K) - 2 * cl.sum_p0(K) / K
                         + cl.N() / (K * K))
    Sfx = K * Sfx
    # gamma(x) = gA - gB x
    one_m_q = rf(1) - q
    gA = c / one_m_q
    gB = 1 / one_m_q
    # Phi(x) = K(gA - gB x) GE - E1 = d x + f  with d = -K gB GE (the
    # x-coefficient, positive here), f = K gA GE - E1
    d = (K * gB * GE) * (-1)
    f_ = K * gA * GE - E1
    # Q(x) = A x^2 + B x + C
    QA = K * (gB * gB * G2) + K
    QB = K * (-2 * gA * gB * G2) + 2 * gB * G1 - 2
    QC = (K * gA * gA * G2 - 2 * gA * G1 + K * (m - 1) / (K * K)
          + K / (K * K) + Sfx)
    return {"K": K, "m": m, "c": c, "q": q, "u0": u0, "r": r,
            "U0": U0, "U1": U1, "U2": U2, "x2": x2, "E2": E2,
            "QA": QA, "QB": QB, "QC": QC, "d": d, "f": f_}


def claims_for(D):
    """Square-completion claims (review-10 route), all (RF, strict):
      C0a/C0b: prefix U(0), U(1) >= U(2);
      S1: A > 0; S2: H2 = E - d^2/A > 0;
      S3: -H'(3/4) > 0 with H' = 2 H2 z + H1, H1 = 2f - Bd/A;
      S4: H(3/4) >= 0.
    """
    A, B, C = D["QA"], D["QB"], D["QC"]
    d, f_ = D["d"], D["f"]
    E = D["K"] * D["E2"]
    H2 = E - d * d / A
    H1 = 2 * f_ - B * d / A
    Hc = C - D["U2"] - B * B / (4 * A)
    z34 = rf_frac(F(3, 4))
    claims = {
        "C0a_U0-U2>=0": (D["U0"] - D["U2"], False),
        "C0b_U1-U2>=0": (D["U1"] - D["U2"], False),
        "S1_A>0": (A, True),
        "S2_H2>0": (H2, True),
        "S3_-Hp(3/4)>0": ((H1 + 2 * H2 * z34) * (-1), True),
        "S4_H(3/4)>=0": (H2 * (z34 * z34) + H1 * z34 + Hc, False),
    }
    return claims


# ---------------- exact per-K reference ----------------
def numeric_reference(K):
    """Exact per-K values from replay: prefix U's, x2/x3, identity check
    values for the quadratic anchors."""
    P0 = replay.exact_bin_probs(replay.near_cdf_exact, K)
    g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
    m = replay.anchor_m_b(K)[0]
    states = [P0]
    sels = []
    for t in range(8):
        D = replay.cumsum_excess(states[-1], K)
        s = replay.argmax_first(D)
        sels.append(s)
        states.append(replay.update_state(states[-1], g, s, ALPHA))
    U = [K * sum((y - F(1, K)) ** 2 for y in states[t]) for t in range(4)]
    return {"sels": sels, "m": m, "U": U,
            "x2": states[2][m - 1], "x3": states[3][m - 1],
            "s1": sels[1], "t2": sels[2],
            "C0": U[0] >= U[2] and U[1] >= U[2]}


def verify_prefix_pattern(nmax=N_CHECK):
    """Replay the canonical CYCLE prefix (s0, s1-closed, t*=2,
    residence) for every CYCLE K in [9, 6*nmax] except the finite
    difference set (covered by the D certificates)."""
    for K in range(9, 6 * nmax + 1):
        if K % 6 not in (3, 4, 5):
            continue
        if K in FINITE_EXCEPTIONS:
            continue
        P = replay.exact_bin_probs(replay.near_cdf_exact, K)
        g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
        m = replay.anchor_m_b(K)[0]
        p = P
        for t in range(6):
            D = replay.cumsum_excess(p, K)
            s = replay.argmax_first(D)
            want = (K // 2 if t == 0 else
                    s1_closed(K) if t == 1 else None)
            if t == 2 and s not in (m - 1, m):
                return (K, t, s, "boundary")
            if t >= 2 and s not in (m - 1, m):
                return (K, t, s, "residence")
            if want is not None and s != want:
                return (K, t, s, want)
            p = replay.update_state(p, g, s, ALPHA)
    return None


def lemma_t2_cycle():
    """CYCLE version of the K>=33 entry lemma: same chain, anchor
    D_{2,m-1}.  Review-10 correction: the odd-K M_0 lower bound is
    M(33) = 200/1089 (K = 33 and 35 are CYCLE; the LOCK lemma's M(37)
    does not carry over).  f1(M33) = 255907/1581228 < 1/6."""
    fails = []
    sixth = F(1, 6)
    M33 = F(3, 16) - F(1, 8 * 33) - F(1, 16 * 33 * 33)
    f1 = lambda M: (F(1, 2) + 3 * M / 4) ** 2 - 4 * M / 3
    f2 = lambda M: (F(1, 2) + M) ** 2 - 5 * M / 3
    d1 = lambda M: F(9, 8) * M - F(7, 12)
    d2 = lambda M: 2 * M - F(2, 3)
    if not M33 == F(200, 1089):
        fails.append("M(33) != 200/1089")
    if not (d1(M33) < 0 and d1(F(3, 16)) < 0 and d2(M33) < 0
            and d2(F(3, 16)) < 0):
        fails.append("f1/f2 not decreasing")
    if not (f1(M33) == F(255907, 1581228) and f1(M33) < sixth
            and f1(F(3, 16)) == F(657, 4096)):
        fails.append("f1 bound")
    if not (f2(M33) < sixth and f2(F(3, 16)) == F(656, 4096)):
        fails.append("f2 bound")
    if not (F(11, 96) < F(9, 64) < sixth and F(5, 36) < F(9, 64)):
        fails.append("non-peak uniform bound")
    if not (F(41, 64) + F(31, 64 * 33) < F(2, 3)):
        fails.append("s_1/K >= 2/3")
    # anchor: D_{2,m-1} > 1/8 for K >= 33 (band-3 top bin, x >= 5/6-1/33,
    # D2 decreasing in Sigma on band 3 so lower-bound with Sigma <= 7/16)
    S_hi = F(7, 16)
    x_lo = F(5, 6) - F(1, 33)
    D2lb = -x_lo * x_lo + (1 + 2 * S_hi) * x_lo - 5 * S_hi / 3
    if not (D2lb > F(1, 8)):
        fails.append(f"anchor D_(2,m-1) <= 1/8 ({float(D2lb)})")
    M1lb = (F(55, 64) - F(3, 16)) / 4
    if not (F(5, 8) <= F(41, 64) - F(33, 64 * 33)):
        fails.append("s1/K lower bound < 5/8")
    if not (M33 + M1lb > F(1, 3)):
        fails.append("M_0 + M_1 <= 1/3")
    if not (F(1, 2) + 3 * F(1, 3) / 4 == F(3, 4) and F(3, 4) > F(2, 3)
            and F(1, 2) + F(1, 3) >= F(5, 6)):
        fails.append("vertex bounds")
    if not (1 - 2 * F(5, 6) < 0):
        fails.append("suffix monotonicity")
    return fails


# ---------------- certificates ----------------
def certify(rf_, n0, strict):
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


def verify_finite_exceptions_via_d():
    """The finite difference set {9,10,11,16,17,23} is covered by the
    COMMITTED D certificates (all-time).  Re-verify exactly those records
    through the INDEPENDENT D verifier (second implementation) with its
    certificate path patched to the subset.  Returns failure list."""
    import importlib.util
    import io
    import contextlib
    import tempfile
    fails = []
    dpath = (REPO / "experiments" / "v03_d_stopping" / "results" /
             "d_certificates.json")
    records = json.loads(dpath.read_text(encoding="utf-8"))
    subset = [r for r in records if r["K"] in FINITE_EXCEPTIONS]
    if sorted(r["K"] for r in subset) != FINITE_EXCEPTIONS:
        return [f"d_certificates.json does not contain exactly the "
                f"finite set {FINITE_EXCEPTIONS}"]
    spec = importlib.util.spec_from_file_location(
        "vdi_e", str(REPO / "experiments" / "v03_d_stopping" /
                     "verify_d_independent.py"))
    vdi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vdi)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "subset.json"
        tmp.write_text(json.dumps(subset), encoding="utf-8")
        old = vdi.CERTS
        old_exp = vdi.expected_ks
        vdi.CERTS = type(old)(tmp)
        vdi.expected_ks = lambda: set(FINITE_EXCEPTIONS)
        try:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = vdi.main()
        finally:
            vdi.CERTS = old
            vdi.expected_ks = old_exp
    if rc != 0:
        fails.append("independent D verifier FAILED on the finite-set "
                     f"subset: {buf.getvalue()[-300:]}")
    return fails


def main():
    t0 = time.time()
    print("== E: CYCLE all-K all-time stopping certificates (review-10 "
          "repaired: 2zPhi identity + square completion) ==")
    lt2 = lemma_t2_cycle()
    if lt2:
        raise SystemExit(f"[E] ABORT: lemma_t2_cycle failed: {lt2}")
    print("[E] lemma t*=2 (CYCLE, K>=33, anchor D_{2,m-1}, M(33) bound) "
          "verified")
    bad = verify_prefix_pattern()
    if bad is not None:
        raise SystemExit(f"[E] ABORT: CYCLE prefix violated at K={bad[0]}"
                         f" step {bad[1]}: {bad}")
    print(f"[E] CYCLE prefix (s0, s1-closed, t*=2, residence) "
          f"replay-verified for every CYCLE K in [9, {6 * N_CHECK}] "
          f"except the finite set {FINITE_EXCEPTIONS}")
    dfe = verify_finite_exceptions_via_d()
    if dfe:
        raise SystemExit(f"[E] ABORT: finite-exception D verification "
                         f"failed: {dfe}")
    print(f"[E] finite set {FINITE_EXCEPTIONS} re-verified through the "
          "INDEPENDENT D verifier (all-time certificates)")
    out = {"families": [],
           "finite_exceptions": [
               {"K": k, "via": "d_certificates.json subset re-verified by "
                "verify_d_independent (all-time: prefix + license + "
                "infinite tail)"} for k in FINITE_EXCEPTIONS],
           "config": {"experiment": "v03_e_cycle_allk",
                      "identity": "U = Q + 2 z Phi + K E2 z^2",
                      "route": "square completion (S1-S4)",
                      "prefix_check_K_max": 6 * N_CHECK}}
    subs = subfamilies()
    for sub in subs:
        r, nu, u0 = sub["r"], sub["nu"], sub["u0"]
        D = build_family(sub)
        # anchors: U-values AND the corrected identity by direct
        # back-substitution at t = 2 (z = 1) and t = 3 (z = 3/4)
        for u_ in (u0, u0 + 1, u0 + 7):
            n = 32 * u_ + nu
            K = 6 * n + r
            ref = numeric_reference(K)
            if ref["s1"] != s1_closed(K) or ref["t2"] not in \
                    (ref["m"] - 1, ref["m"]):
                raise SystemExit(f"[E] ABORT r={r} nu={nu} K={K}: "
                                 f"prefix mismatch {ref['s1']}")
            if not ref["C0"]:
                raise SystemExit(f"[E] ABORT r={r} nu={nu} K={K}: C0")
            for key in ("U0", "U1", "U2"):
                symv = D[key]
                val = pp_eval(symv.num, u_) / pp_eval(symv.den, u_)
                if val != ref["U"][int(key[1])]:
                    raise SystemExit(f"[E] ABORT r={r} nu={nu} K={K}: "
                                     f"{key} symbolic {val} != "
                                     f"{ref['U'][int(key[1])]}")
            # DIRECT BACK-SUBSTITUTION (review-10 mandate): U(p_t) equals
            # the assembled quadratic Q(x_t) + 2 z Phi(x_t) + K E2 z^2
            A_ = pp_eval(D["QA"].num, u_) / pp_eval(D["QA"].den, u_)
            B_ = pp_eval(D["QB"].num, u_) / pp_eval(D["QB"].den, u_)
            C_ = pp_eval(D["QC"].num, u_) / pp_eval(D["QC"].den, u_)
            d_ = pp_eval(D["d"].num, u_) / pp_eval(D["d"].den, u_)
            f_ = pp_eval(D["f"].num, u_) / pp_eval(D["f"].den, u_)
            E2_ = pp_eval(D["E2"].num, u_) / pp_eval(D["E2"].den, u_)
            for t, xt in ((2, ref["x2"]), (3, ref["x3"])):
                z = RHO ** (t - 2)
                quad = (A_ * xt * xt + B_ * xt + C_) + \
                    2 * z * (d_ * xt + f_) + K * E2_ * z * z
                if quad != ref["U"][t]:
                    raise SystemExit(
                        f"[E] ABORT r={r} nu={nu} K={K}: identity "
                        f"back-substitution FAILS at t={t} "
                        f"(quad {quad} != U {ref['U'][t]})")
        certs = {}
        for name, (rf_, strict) in claims_for(D).items():
            st = certify(rf_, u0, strict)
            if st is None:
                raise SystemExit(f"[E] ABORT r={r} nu={nu} {name}: no "
                                 f"certificate (num degree "
                                 f"{len(rf_.num) - 1})")
            vf = verify_claim(rf_, st, u0, strict)
            if vf:
                raise SystemExit(f"[E] ABORT r={r} nu={nu} {name}: {vf}")
            certs[name] = {"num": st["num"], "den_sign": st["den_sign"],
                           "den_cert": st["den_cert"]}
        out["families"].append({
            "r": r, "nu": nu, "u0": u0, "claims": certs,
            "num_degrees": {k: len(p.num) - 1
                            for k, (p, _) in claims_for(D).items()},
        })
        if len(out["families"]) % 16 == 0:
            print(f"[E] {len(out['families'])}/96 sub-families certified")
    out["git"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                                capture_output=True, text=True).stdout.strip()
    out["source_sha256"] = _source_sha256()
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "e_cycle_cert.json").write_text(json.dumps(out, indent=1),
                                               encoding="utf-8")
    fails = verify_from_json()
    if fails:
        raise SystemExit(f"[E] ABORT: JSON re-verification failed: "
                         f"{fails[:5]}")
    print(f"[E] PASS: 96/96 sub-families (square-completion route) + "
          f"finite set via D certificates + JSON re-verified "
          f"({time.time() - t0:.1f}s)")
    return 0


def verify_from_json():
    data = json.loads((RESULTS / "e_cycle_cert.json").read_text(
        encoding="utf-8"))
    fails = []
    if data.get("source_sha256") != _source_sha256():
        fails.append("source_sha256 mismatch")
    if data.get("config", {}).get("identity") != "U = Q + 2 z Phi + K E2 z^2" \
            or data.get("config", {}).get("route") != \
            "square completion (S1-S4)":
        fails.append("config does not record the review-10 identity/route")
    expected = {(s["r"], s["nu"], s["u0"]) for s in subfamilies()}
    got = [(f["r"], f["nu"], f["u0"]) for f in data["families"]]
    if sorted(got) != sorted(expected):
        fails.append("sub-family set mismatch (or start ranges)")
    for entry in data["families"]:
        key = (entry["r"], entry["nu"])
        sub = next((s for s in subfamilies()
                    if (s["r"], s["nu"]) == key), None)
        if sub is None or entry["u0"] != sub["u0"]:
            continue
        D = build_family(sub)
        claims = claims_for(D)
        stored = set(entry["claims"])
        if stored != set(claims):
            fails.append(f"r={entry['r']} nu={entry['nu']}: claim set "
                         f"mismatch")
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
    if [e["K"] for e in data.get("finite_exceptions", [])] != \
            FINITE_EXCEPTIONS:
        fails.append(f"finite exception set != {FINITE_EXCEPTIONS}")
    elif not all("verify_d_independent" in e.get("via", "") and
                 "d_certificates.json" in e.get("via", "")
                 for e in data["finite_exceptions"]):
        fails.append("finite exceptions must reference the D certificates "
                     "re-verified through verify_d_independent")
    else:
        fails += verify_finite_exceptions_via_d()
    return fails


if __name__ == "__main__":
    raise SystemExit(main())
