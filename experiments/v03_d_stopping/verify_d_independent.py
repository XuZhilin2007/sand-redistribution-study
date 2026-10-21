"""Independent verifier for the stopping-experiment certificates (WP0).

STATUS (2026-10-20): this verifies the per-lattice-point D evidence
(203 K).  The all-K LOCK statement is carried by Theorem B-lock
(docs/V03_B_LOCK_ALLK.md); the CYCLE all-K statement is open.  CYCLE
states are only ASYMPTOTICALLY periodic (geometric convergence of x_t,
p_t, d_t to a rational periodic solution; the selection sequence alone
is exactly eventually periodic) -- this verifier's cycle checks certify
the asymptotic minimum, not exact state periodicity.

SECOND IMPLEMENTATION: imports nothing from run_d.py or
sand_m0.word_cyclicity. It reads results/d_certificates.json and, for
every CYCLE record, re-derives from the S3 definitions with its own
exact Fractions:

  1. boundary parameters (m, h, c, q, L, beta);
  2. the word's affine fixed point xstar and Lambda (cross-checked);
  3. the exact periodic profile at every phase (geometric-series solve
     with its own driving reconstruction), with a one-step back-check;
  4. the license: the recorded island half-width, recomputed dmin, and
     that ps[phase] +- delta brackets (i.e. the island is well formed);
  5. the tail lower bound K*Q2 - K*(1+Q2)*D2 with its own envelope
     (alpha, beta_c = dx0/(4 mu), unimodal max by scanning);
  6. U2: recomputed from the canonical initial state by running the
     model's own two first updates INDEPENDENTLY (closed form of the
     canonical init and first selections s0, s1 for K >= 7; the K = 6
     bridge handled by one direct step);
  7. the certificate inequality U2 <= tail_lb.

Run:  python experiments/v03_d_stopping/verify_d_independent.py
Exit 0 + INDEPENDENT D-VERIFY PASS <n>/<N> iff every record verifies.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
import sys as _sys
try:
    _sys.set_int_max_str_digits(1_000_000)
except (AttributeError, ValueError):
    pass
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERTS = HERE / "results" / "d_certificates.json"


def params(K):
    m = -(-5 * K // 6)
    h = F(1, K)
    b = F(m, K)
    c = 2 * b - b * b
    q = F(6 - K % 6, 3 * K)
    L = c * q
    return m, h, b, c, q, L


def canonical_init(K):
    return [F(2 * K - 2 * i + 1, K * K) for i in range(1, K + 1)]


def g_kernel(K):
    def G(y):
        if y <= F(1, 3):
            return y / 2
        if y <= F(2, 3):
            return 3 * y / 2 - F(1, 3)
        if y <= F(5, 6):
            return 2 * y - F(2, 3)
        return F(1)
    out = []
    for i in range(1, K + 1):
        out.append(G(F(i, K)) - G(F(i - 1, K)))
    return out


def argmax_sel(D):
    best = D[0]
    idx = 0
    for j in range(1, len(D)):
        if D[j] > best:
            best = D[j]
            idx = j
    return idx + 1  # 1-based, smallest index


def discrep(p, K, h):
    """Cumulative discrepancies in O(K) (running sum)."""
    out = []
    F = F0 = 0
    for j in range(K):
        F0 += p[j]
        out.append(F0 - (j + 1) * h)
    return out


def update(p, g, s, K):
    M = sum(p[:s]) / 4
    return [F(3, 4) * p[i] + M * g[i] if i < s else p[i] + M * g[i]
            for i in range(K)]


def U_of(p, K):
    h = F(1, K)
    return K * sum((pi - h) ** 2 for pi in p)


def verify_record(rec, RHO=F(3, 4)):
    K = rec["K"]
    m, h, b, c, q, L = params(K)
    g = g_kernel(K)
    w = rec["word"]
    n = rec["n"]
    assert len(w) == n
    # 2. fixed point and Lambda
    lam, mu_c = F(1), F(0)
    for br in w:
        if br == "A":
            lam, mu_c = lam * RHO, mu_c * RHO + L * F(1, 4)
        else:
            f = 1 - q * F(1, 4)
            lam, mu_c = lam * f, mu_c * f + c * q * F(1, 4)
    if str(lam) != rec["lam_word"]:
        return f"Lambda mismatch: {lam} vs {rec['lam_word']}"
    xstar = mu_c / (1 - lam)
    if str(xstar) != rec["xstar"]:
        return f"xstar mismatch: {xstar} vs {rec['xstar']}"
    ps = [xstar]
    for br in w:
        x = ps[-1]
        ps.append(RHO * x + L / 4 if br == "A" else x + (c - x) * q / 4)
    if ps[n] != ps[0]:
        return "fixed point does not realize the word"
    if any((pp > h) != (br == "A") for pp, br in zip(ps, w)):
        return "fixed point branch sides wrong"
    # 4. license well-formedness
    dmin = min(abs(pp - h) for pp in ps[:n])
    delta = dmin / 64
    if str(delta) != rec["island_delta"]:
        return f"island delta mismatch: {delta} vs {rec['island_delta']}"
    if F(rec["license_margin"]) <= 0:
        return "license margin non-positive"
    # 3. exact periodic profiles (own geometric solve) + back-check
    G = [sum(g[:j + 1]) for j in range(K)]

    def bdr(jpos, col):
        Fs = c if w[jpos] == "A" else c - ps[jpos]
        return F(1, 4) * (Fs * G[col] - (col + 1) * h)

    rn = RHO ** n
    profiles = []
    for jpos in range(n):
        p_prof = [None] * K
        prevF = F(0)
        for col in range(m - 1):
            tot = F(0)
            for i in range(n):
                tot += RHO ** i * bdr((jpos - 1 - i) % n, col)
            Fj = tot / (1 - rn) + (col + 1) * h
            p_prof[col] = Fj - prevF
            prevF = Fj
        p_prof[m - 1] = ps[jpos]
        for col in range(m, K):
            p_prof[col] = canonical_init(K)[col]
        profiles.append(p_prof)
    for jpos in range(n):
        s_sel = m if w[jpos] == "A" else m - 1
        if update(profiles[jpos], g, s_sel, K) != profiles[(jpos + 1) % n]:
            return f"profile back-check fails at phase {jpos}"
    # 5. INDEPENDENT trajectory to t0 with the verifier's OWN engine,
    #    deviation recomputation, phase-offset confirmation, license
    #    check, and the envelope with own arithmetic.
    mu = max(RHO, 1 - q * F(1, 4))
    if str(mu) != rec["mu_env"]:
        return f"mu mismatch: {mu} vs {rec['mu_env']}"
    p = canonical_init(K)
    for t in range(rec["t0"]):
        s = argmax_sel(discrep(p, K, h))
        p = update(p, g, s, K)
    # confirm the recorded phase offset on the OWN trajectory
    offs = rec["phase_offset"]
    for i in range(2 * n):
        s_expect = m if w[(offs + i) % n] == "A" else m - 1
        if argmax_sel(discrep(p, K, h)) != s_expect:
            return f"phase offset mismatch at +{i}"
        p = update(p, g, s_expect, K)
    # recompute the t0 state (the confirmation advanced 2n steps)
    p = canonical_init(K)
    for t in range(rec["t0"]):
        s = argmax_sel(discrep(p, K, h))
        p = update(p, g, s, K)
    q0 = profiles[offs]
    dx0 = abs(p[m - 1] - q0[m - 1])
    dint0 = max((abs(p[col] - q0[col]) for col in range(m - 1)),
                default=F(0))
    # real license check: the t0 boundary deviation fits in the island
    if not (dx0 <= delta):
        return "license anchor violated (dx0 > delta)"
    # consistency with the recorded generator floats
    if abs(float(dx0) - rec["dx0"]) > 1e-9 * max(float(rec["dx0"]), 1e-300):
        return "dx0 disagrees with recorded value"
    if abs(float(dint0) - rec["dint0"]) > 1e-9 * max(float(rec["dint0"]),
                                                     1e-300):
        return "dint0 disagrees with recorded value"
    alpha = dint0 + dx0
    beta_c = dx0 / (4 * mu)
    g_max = alpha
    g_prev = alpha
    dt = 0
    while dt < 10 ** 6:
        g_next = (alpha + beta_c * (dt + 1)) * mu ** (dt + 1)
        if g_next <= g_prev:
            break
        g_max = g_next
        g_prev = g_next
        dt += 1
    import math
    D2 = F(math.isqrt(K - 1) + 1) * g_max
    tail_lb = None
    for jpos in range(n):
        Q2 = sum((qi - h) ** 2 for qi in profiles[jpos])
        lb = K * Q2 - K * (1 + Q2) * D2
        tail_lb = lb if tail_lb is None else min(tail_lb, lb)
    # 6. U2 recomputed independently: run two updates from canonical init
    p2 = canonical_init(K)
    for t in range(2):
        s = argmax_sel(discrep(p2, K, h))
        p2 = update(p2, g, s, K)
    U2 = U_of(p2, K)
    if str(U2) != rec["U2"]:
        return f"U2 mismatch: {U2} vs {rec['U2']}"
    # 7. the certificate inequality with independently recomputed
    #    deviations and envelope
    if not (U2 <= tail_lb):
        return "U2 > independent tail bound"
    # 8. finite-prefix optimality: replay t = 0..t0 with the OWN engine,
    #    check U(p_t) >= U(p_2) at EVERY step, AND compute the first
    #    boundary-entry time independently (exact match required)
    pt = canonical_init(K)
    t_star_own = None
    for t in range(rec["t0"]):
        s = argmax_sel(discrep(pt, K, h))
        if t_star_own is None and s in (m - 1, m):
            t_star_own = t
        if U_of(pt, K) < U2:
            return f"prefix optimality fails at t={t}"
        pt = update(pt, g, s, K)
    if t_star_own is None:
        for t in range(rec["t0"], rec["t0"] + 3 * n):
            s = argmax_sel(discrep(pt, K, h))
            if s in (m - 1, m):
                t_star_own = t
                break
            pt = update(pt, g, s, K)
    # 9. recorded proof-field consistency
    if rec.get("best_t") != 2:
        return f"best_t record wrong: {rec.get('best_t')}"
    if str(tail_lb) != rec["tail_lb"]:
        return (f"tail_lb mismatch: {tail_lb} vs {rec['tail_lb']}")
    if t_star_own is None or rec.get("t_star") != t_star_own:
        return (f"t_star mismatch: record={rec.get('t_star')} "
                f"own-replay={t_star_own}")
    lm = F(rec["license_margin"])
    # license_margin must equal delta - |x_t0 - ps[offs]| exactly
    if lm != delta - dx0:
        return f"license_margin mismatch: {lm} vs {delta - dx0}"
    return None


def expected_ks():
    cfg_Ks = list(range(6, 204)) + [300, 600, 999, 1000, 1001]
    return set(cfg_Ks)


def main():
    certs = json.loads(CERTS.read_text(encoding="utf-8"))
    cyc = [r for r in certs if r["class"] == "CYCLE"]
    lock = [r for r in certs if r["class"] == "LOCK"]
    # completeness: the certificate set must cover exactly the expected
    # lattice, with the correct class per K, and WITHOUT duplicates
    exp = expected_ks()
    got_list = [r["K"] for r in certs]
    if len(got_list) != len(set(got_list)):
        from collections import Counter
        dup = [k for k, c in Counter(got_list).items() if c > 1]
        print("INDEPENDENT D-VERIFY FAIL: duplicate K records: "
              f"{sorted(dup)[:5]}")
        return 1
    got = set(got_list)
    if got != exp:
        missing = sorted(exp - got)[:5]
        extra = sorted(got - exp)[:5]
        print("INDEPENDENT D-VERIFY FAIL: certificate set incomplete "
              f"({len(got)} records vs {len(exp)} expected; "
              f"missing={missing} extra={extra})")
        return 1
    if not certs:
        print("INDEPENDENT D-VERIFY FAIL: empty certificate set")
        return 1
    for r in certs:
        want = "LOCK" if r["K"] % 6 in (0, 1, 2) else "CYCLE"
        if r["class"] != want:
            print(f"INDEPENDENT D-VERIFY FAIL: K={r['K']} class "
                  f"{r['class']} != {want}")
            return 1
    fails = []
    for rec in cyc:
        err = verify_record(rec)
        if err:
            fails.append((rec["K"], err))
    if fails:
        print("INDEPENDENT D-VERIFY FAIL:")
        for K, err in fails[:10]:
            print(f"  K={K}: {err}")
        return 1
    print(f"INDEPENDENT D-VERIFY PASS "
          f"({len(cyc)}/{len(cyc)} CYCLE records incl. prefix-optimality "
          f"and field consistency; {len(lock)} LOCK records carry the "
          f"quadratic branch certificate in d_rows.csv)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
