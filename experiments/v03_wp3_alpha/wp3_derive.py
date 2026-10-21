"""WP3 — general-alpha LN dictionary: exact derivation & verification.

A  conjugacy + dictionary + eta alpha-independence (exact Fractions;
   pointwise identity checks as in V03_LITERATURE_COMPARISON 0.5.2 —
   both sides affine in x and ratio-of-linear in alpha, so >=2 exact
   points per branch prove the identity; we use 26 points x 6 alphas).
B  LN standing-hypothesis region on the (K, alpha) grid K=9..33.
C  existence boundary delta = 1-lambda  <=>  eta = 1  <=>  L = h.
D  plateau vs exact-itinerary spot checks (entry via full vector,
   then exact scalar map; pattern copied from run_c.entry_and_itinerary).
"""
from __future__ import annotations

import sys
from fractions import Fraction as F
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "s4_w1_s1_xo_precision"))
import replay  # noqa: E402

ALPHAS_ETA = [F(1, 4), F(1, 8), F(1, 2), F(3, 4), F(2, 3), F(1, 100)]
KS_ETA = [9, 10, 11, 21, 65]
ALPHAS_GRID = [F(1, 100), F(1, 8), F(1, 6), F(1, 4), F(1, 3),
               F(1, 2), F(2, 3), F(3, 4), F(9, 10)]
Q_MAX = 64


def dictionary(K, alpha):
    """Exact general-alpha dictionary on U=(u_-,u_+]."""
    m = (5 * K + 5) // 6
    h = F(1, K)
    b = F(m, K)
    c = 2 * b - b * b
    q = F(6 - K % 6, 3 * K)
    L = c * q
    A = lambda x: (1 - alpha) * x + alpha * L          # x > h
    B = lambda x: x + alpha * (c - x) * q               # x <= h
    um, up = A(h), B(h)
    D = up - um
    lam = 1 - alpha
    mu = (1 - alpha * q) / (1 - alpha)
    eta = (up - h) / D
    delta = 1 - lam * eta
    return {"K": K, "m": m, "h": h, "b": b, "c": c, "q": q, "L": L,
            "u_minus": um, "u_plus": up, "D": D,
            "lam": lam, "mu": mu, "eta": eta, "delta": delta}


# ---------------------------------------------------------------- Part A
def part_a():
    print("=== A: conjugacy + dictionary + eta(alpha)-independence ===")
    for K in KS_ETA:
        etas = set()
        for alpha in ALPHAS_ETA:
            d = dictionary(K, alpha)
            h, q, c, L = d["h"], d["q"], d["c"], d["L"]
            um, up, D = d["u_minus"], d["u_plus"], d["D"]
            assert D == alpha * h * (1 - q), (K, alpha, "D=alpha h(1-q)")
            assert up == (1 - alpha * q) * h + alpha * L, (K, alpha, "u+")
            assert um == (1 - alpha) * h + alpha * L, (K, alpha, "u-")
            assert um < h < up, (K, alpha, "U straddles h")
            # raw conjugacy constant must coincide with delta=1-lambda*eta
            assert alpha * (up - L) / D == d["delta"], (K, alpha, "delta raw")
            assert d["eta"] == q * (c - h) / ((1 - q) * h), (K, alpha, "eta")
            etas.add(d["eta"])
            # pointwise branch identity at 26 exact rational points
            pts = [um + F(i, 25) * (up - um) for i in range(1, 26)]
            pts += [h, up]
            for x in pts:
                assert um < x <= up
                z = (up - x) / D
                if x > h:
                    xp = (1 - alpha) * x + alpha * L
                    want = d["lam"] * z + d["delta"]
                else:
                    xp = x + alpha * (c - x) * q
                    want = d["mu"] * (d["lam"] * z + d["delta"] - 1)
                zp = (up - xp) / D
                assert zp == want, (K, alpha, x, "branch identity")
                assert um < xp <= up, (K, alpha, x, "U invariant")
        assert len(etas) == 1, (K, "eta depends on alpha!")
        eta = etas.pop()
        d0 = dictionary(K, F(1, 4))
        print(f"K={K:3d}: eta={eta} identical for all 6 alphas; "
              f"branch identity + U-invariance OK at 28 pts x 6 alphas")
    print("A: CONFIRMED lambda=1-alpha, mu=(1-alpha q)/(1-alpha), "
          "delta=1-(1-alpha)eta, eta alpha-free")


# ---------------------------------------------------------------- Part B
def part_b():
    print("\n=== B: LN standing hypotheses on the grid ===")
    print("     (0<lam<1, mu>0, lam*mu<1 => d=1, 1-lam<delta<1)")
    labels = ["1/100", "1/8", "1/6", "1/4", "1/3", "1/2", "2/3", "3/4", "9/10"]
    print("K    r class  eta                " + " ".join(f"{s:>5}" for s in labels))
    for K in range(9, 34):
        row = []
        eta = None
        for alpha in ALPHAS_GRID:
            d = dictionary(K, alpha)
            hyp = (0 < d["lam"] < 1 and d["mu"] > 0
                   and d["lam"] * d["mu"] < 1
                   and 1 - d["lam"] < d["delta"] < 1)
            row.append(hyp)
            # exact margin identities
            assert d["lam"] * d["mu"] == 1 - alpha * d["q"]
            assert d["delta"] - (1 - d["lam"]) == (1 - alpha) * (1 - d["eta"])
            assert d["delta"] < 1 or d["eta"] <= 0
            eta = d["eta"]
        assert len(set(row)) == 1, (K, "boolean varies with alpha!")
        cyc = K % 6 in (3, 4, 5)
        assert row[0] == cyc, (K, "boolean != CYCLE class")
        Lh = dictionary(K, F(1, 4))["L"] / dictionary(K, F(1, 4))["h"]
        assert Lh < 1 if cyc else Lh > 1
        print(f"{K:3d}  {K % 6} {'CYCLE' if cyc else 'LOCK '}  "
              f"{str(eta):18s} " + " ".join(f"{'Y' if v else 'n':>5}" for v in row)
              + f"   L/h={Lh}")
    print("B: hypotheses hold exactly on CYCLE classes (r in {3,4,5}) "
          "for EVERY alpha in the grid; never on LOCK; alpha-independent.")


# ---------------------------------------------------------------- Part C
def part_c():
    print("\n=== C: existence / degeneracy boundary ===")
    hits = []
    for K in range(6, 2001):
        d = dictionary(K, F(1, 4))
        if d["eta"] == 1 or d["L"] == d["h"]:
            hits.append(K)
        assert (d["eta"] == 1) == (d["L"] == d["h"])
        assert (d["eta"] < 1) == (d["L"] < d["h"])
        assert (d["u_minus"] < d["h"]) == (d["L"] < d["h"])
    assert not hits, hits
    print("C: no K in 6..2000 has eta=1 (delta=1-lambda); U nonempty "
          "iff L<h iff c(6-r)/3<1 iff r in {3,4,5} — alpha-independent.")


# ---------------------------------------------------------------- Part D
def S_sum(lam, mu, p, q):
    if q == 2:
        return F(0)
    tot = F(0)
    for k in range(1, q - 1):
        dd = ((k + 1) * p) // q - (k * p) // q
        tot += dd * lam ** k * mu ** ((k * p) // q)
    return tot


def delta_plus(lam, mu, p, q):
    S = S_sum(lam, mu, p, q)
    num = (1 - lam) * (1 + mu * S + lam ** (q - 1) * mu ** p * (1 - lam))
    den = 1 + (mu - 1) * S + lam ** (q - 1) * mu ** (p - 1) * \
        (mu - lam * mu - 1)
    return num / den


def delta_minus(lam, mu, p, q):
    S = S_sum(lam, mu, p, q)
    num = (1 - lam) * (1 + mu * S)
    den = 1 + (mu - 1) * S - lam ** q * mu ** (p - 1)
    return num / den


def find_plateau(lam, mu, delta):
    for q in range(2, Q_MAX + 1):
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            dlo = delta_minus(lam, mu, p, q)
            dhi = delta_plus(lam, mu, p, q)
            if dlo <= delta <= dhi:
                return p, q, dlo, dhi
    return None


def entry_and_itinerary(K, alpha, H, full_head=150):
    """Exact full vector until entry+full_head (scalar closure asserted),
    then exact scalar map only. Returns (t_star, selection list t=0..H-1)."""
    m = (5 * K + 5) // 6
    h = F(1, K)
    b = F(m, K)
    c = 2 * b - b * b
    q = F(6 - K % 6, 3 * K)
    L = c * q
    p = replay.exact_bin_probs(replay.near_cdf_exact, K)
    g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
    s_list = []
    t_star = None
    x = None
    t = 0
    while t < H:
        if t_star is None or t <= t_star + full_head:
            D = replay.cumsum_excess(p, K)
            s = replay.argmax_first(D)
            if t_star is not None:
                assert p[m - 1] == x, f"scalar closure failed K={K} t={t}"
            if t_star is None and s in (m - 1, m):
                t_star = t
                x = p[m - 1]
            if t_star is not None:
                x = ((1 - alpha) * x + alpha * L if x > h
                     else x + alpha * (c - x) * q)
            p = replay.update_state(p, g, s, alpha)
        else:
            s = m if x > h else m - 1
            x = ((1 - alpha) * x + alpha * L if x > h
                 else x + alpha * (c - x) * q)
        s_list.append(s)
        t += 1
    return t_star, s_list


def eventual_period_word(letters):
    """letters: 'A'/'B' string, long enough that the last 1500 are purely
    periodic. Returns (period-minimal word, minimal period)."""
    w = letters[-1500:]
    n = len(w)
    for T in range(1, n):
        if w[:n - T] == w[T:]:
            return w[:T], T
    raise RuntimeError("no period in tail")


def mechanical_word(p, q):
    return "".join("A" if (j * p) % q < q - p else "B" for j in range(q))


def same_rotation(w1, w2):
    if len(w1) != len(w2):
        return False
    return any(w1[i:] + w1[:i] == w2 for i in range(len(w1)))


def compress(w):
    out, i = [], 0
    while i < len(w):
        j = i
        while j < len(w) and w[j] == w[i]:
            j += 1
        out.append(f"{w[i]}^{j - i}" if j - i > 1 else w[i])
        i = j
    return "".join(out)


def part_d():
    print("\n=== D: plateau vs exact itinerary (spot checks) ===")
    pairs = [(9, F(1, 4)), (9, F(1, 8)), (9, F(1, 2)), (10, F(1, 8)),
             (10, F(1, 2)), (21, F(1, 8)), (21, F(1, 2))]
    n_checked = 0
    for K, alpha in pairs:
        d = dictionary(K, alpha)
        pl = find_plateau(d["lam"], d["mu"], d["delta"])
        if pl is None:
            raise SystemExit(
                f"WP3 FAIL: K={K} alpha={alpha}: no plateau with "
                f"q<={Q_MAX} — check INCOMPLETE, not passing")
        p, q, dlo, dhi = pl
        T_alpha = int(F(81 * K) / (80 * alpha)) + 1
        H = T_alpha + 150 + 3400
        t_star, s_list = entry_and_itinerary(K, alpha, H)
        assert t_star is not None and t_star <= T_alpha, (K, t_star)
        letters = "".join("B" if s == d["m"] - 1 else "A"
                          for s in s_list[t_star:])
        assert set(letters) == {"A", "B"}, "both branches must appear"
        word, T = eventual_period_word(letters)
        freq = F(word.count("B"), T)
        inter = float((d["delta"] - dlo) / (dhi - dlo))
        assert freq == F(p, q), (K, alpha, freq, p, q)
        assert T == q, (K, alpha, T, q)
        assert dlo < d["delta"] < dhi, "endpoint (iii) excluded"
        mw_ok = same_rotation(word, mechanical_word(p, q))
        assert mw_ok, (K, alpha, "word is NOT a rotation of the "
                                 "mechanical word")
        print(f"K={K:3d} alpha={alpha}: lam={d['lam']} mu={d['mu']} "
              f"delta={d['delta']}")
        print(f"   t*={t_star} word={compress(word)} period={T} "
              f"rho=B-freq={p}/{q}; plateau interior={inter:.4f}; "
              f"mechanical-word rotation match={mw_ok}")
        n_checked += 1
    if n_checked != len(pairs):
        raise SystemExit(f"WP3 INCOMPLETE: {n_checked}/{len(pairs)}")
    print("D: itinerary minimal period == plateau denominator q and "
          "B-frequency == p/q for every spot pair; delta strictly "
          "interior (LN case (ii), endpoint (iii) never hit).")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    part_d()
    print("\nWP3 ALL CHECKS PASSED")
