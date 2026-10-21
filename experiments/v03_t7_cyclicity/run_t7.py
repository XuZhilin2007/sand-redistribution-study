"""v0.3 T7 runner: exact eventual-periodicity of the canonical XO boundary
itinerary (CYCLE classes) and explicit LOCK asymptotics.

Frozen experiment (config.json). All checks are exact Fraction arithmetic.

  Step 0  reduction gate: full-vector exact replay vs the T4 scalar map,
          from first boundary entry, on the census lattice + anchors.
  Step 1  strong word detection (>= REP exact run-period repeats) +
          interval-invariance certificate per lattice K = 9..603,
          r in {3,4,5}, + anchors.
  Step 2  uniform layer:
    2a r=3 pure-window membership per lattice K (global-island criterion),
        sigma monotonicity envelope, anchor windows;
    2b r=4 global cell-alternation inequalities per lattice K in
        [K0_4, K_UNI] + anchor;
    2c r=5 fixed-word fast certificates per lattice K in [K0_5, K_UNI],
        limit-island exact interval covering (BFS), tail shadowing
        constants via exact polynomial expansion.
  Step 3  LOCK lattice sample: last-B time, then s=m, d=a exactly,
          x monotone to L, suffix frozen, exact affine error formula.

Hard halts stop the run; scientific events are recorded, never repaired.
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "s4_w1_s1_xo_precision"))

import replay  # noqa: E402
from sand_m0.word_cyclicity import (  # noqa: E402
    ALPHA, boundary_params, certify_word, compose_word, detect_eventual_word,
    forced_step, pure_window_membership, runs_of, scalar_orbit, scalar_step,
    sigma_strictly_monotone,
)

CFG = json.loads((HERE / "config.json").read_text(encoding="utf-8"))
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

HARD = []
EVENTS = []
K_TAIL = [None]  # set by step2c_tail; finite-K BFS covers [65, K_TAIL)


def halt(reason):
    HARD.append(reason)
    raise SystemExit(f"HARD HALT: {reason}")


def event(kind, detail):
    EVENTS.append({"kind": kind, "detail": detail})


def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_sha():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                              capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:
        return "unknown"


def write_csv(name, fields, rows):
    p = OUT / name
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return p


def lattice(r, lo, hi):
    return [K for K in range(lo, hi + 1) if K % 6 == r]


def census_lattice():
    return (lattice(3, CFG["census_K_lo"], CFG["census_K_hi"])
            + lattice(4, CFG["census_K_lo"], CFG["census_K_hi"])
            + lattice(5, CFG["census_K_lo"], CFG["census_K_hi"])
            + [K for K in CFG["census_anchors"]])


def win3_lattice():
    return lattice(3, CFG["win3_K_lo"], CFG["win3_K_hi"]) + CFG["win3_anchors"]


def alt4_lattice():
    return (lattice(4, CFG["alt4_K_lo"], CFG["alt4_K_hi"])
            + [K for K in CFG["alt4_anchors"]])


def cert5_lattice():
    return lattice(5, CFG["cert5_K_lo"], CFG["cert5_K_hi"])


# ---------------------------------------------------------------------------
# entry state via full exact vector replay (S1/S2 semantics)
# ---------------------------------------------------------------------------

def entry_state_exact(K: int):
    P = boundary_params(K)
    m = P["m"]
    p = replay.exact_bin_probs(replay.near_cdf_exact, K)
    g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
    for t in range(CFG["entry_cap"]):
        D = replay.cumsum_excess(p, K)
        s = replay.argmax_first(D)
        if s in (m - 1, m):
            return t, p[m - 1], p, g
        p = replay.update_state(p, g, s, ALPHA)
    halt(f"K={K}: no boundary entry within entry_cap")


# ---------------------------------------------------------------------------
# Step 0: scalar-reduction gate
# ---------------------------------------------------------------------------

def step0():
    bad = []
    checked = 0
    for K in CFG["step0_K"]:
        if K % 6 not in (3, 4, 5):
            continue
        t_star, x0, p, g = entry_state_exact(K)
        P = boundary_params(K)
        # continue the full vector and the scalar map in parallel
        x = x0
        for t in range(CFG["step0_steps"]):
            D = replay.cumsum_excess(p, K)
            s_full = replay.argmax_first(D)
            x_pred, br_pred = scalar_step(x, P)
            s_pred = P["m"] if x > P["h"] else P["m"] - 1
            if s_full != s_pred or p[P["m"] - 1] != x:
                bad.append((K, t_star + t, s_full, s_pred))
                break
            p = replay.update_state(p, g, s_full, ALPHA)
            x = x_pred
        else:
            checked += 1
    if bad:
        halt(f"step0 reduction mismatch: {bad[:3]}")
    return {"checked_K": checked}


# ---------------------------------------------------------------------------
# Step 1: census + certificates
# ---------------------------------------------------------------------------

def step1():
    rows = []
    for K in census_lattice():
        if K % 6 not in (3, 4, 5):
            continue
        P = boundary_params(K)
        t_star, x0, _, _ = entry_state_exact(K)
        xs, it = scalar_orbit(x0, P, CFG["orbit_T"])
        i0, n, w = detect_eventual_word(it, rep=CFG["rep"])
        if w is None:
            halt(f"K={K}: no >= {CFG['rep']}-repeat word within T="
                 f"{CFG['orbit_T']} (census)")
        cert = certify_word(K, x0, w, xs=xs, it=it, i0=i0,
                            delta_frac=Fraction(*CFG["delta_frac"]))
        if not cert.get("ok"):
            halt(f"K={K}: certificate failed: {cert.get('reason')}")
        rows.append({
            "K": K, "r": K % 6, "m": P["m"], "t_star": t_star,
            "t0": i0, "t_entry": cert["t_entry"], "n": n, "word": w,
            "lam_num": cert["lam"].numerator, "lam_den": cert["lam"].denominator,
            "dmin_num": cert["dmin"].numerator,
            "dmin_den": cert["dmin"].denominator,
            "xstar_num": cert["xstar"].numerator,
            "xstar_den": cert["xstar"].denominator,
            "cycle": "|".join(f"{p.numerator}/{p.denominator}"
                              for p in cert["cycle"]),
            "phase": cert["phase"],
        })
    p = write_csv("t7_census_certificates.csv",
                  list(rows[0].keys()), rows)
    return {"n_rows": len(rows), "csv": p.name,
            "max_t_entry": max(r["t_entry"] for r in rows)}


# ---------------------------------------------------------------------------
# Step 2a: r=3 pure windows
# ---------------------------------------------------------------------------

def step2a():
    rows = []
    for K in win3_lattice():
        if K % 6 != 3:
            continue
        win = pure_window_membership(K)
        rows.append({"K": K, "window_k": win if win else -1})
    mono = sigma_strictly_monotone(CFG["win3_mono_lo"], CFG["win3_mono_hi"], 3)
    if not mono:
        halt("r=3 sigma not strictly monotone on probe lattice")
    p = write_csv("t7_r3_windows.csv", ["K", "window_k"], rows)
    return {"rows": len(rows), "csv": p.name,
            "K_min_pure8": min(r["K"] for r in rows if r["window_k"] == 8),
            "mixed_K": [r["K"] for r in rows if r["window_k"] == -1]}


# ---------------------------------------------------------------------------
# Step 2b: r=4 global cell alternation
# ---------------------------------------------------------------------------

def step2b():
    F34 = Fraction(3, 4), Fraction(4, 3), Fraction(9, 16)
    worst1 = worst2 = None
    n_ok = 0
    for K in alt4_lattice():
        if K % 6 != 4:
            continue
        P = boundary_params(K)
        beta, s = P["beta"], P["sigma"]
        thr = (Fraction(4, 3) - s) / beta
        m1 = Fraction(3, 4) * (beta * Fraction(3, 4) + s) - thr
        m2 = thr - Fraction(9, 16) * (beta + s)
        if not (m1 > 0 and m2 >= 0):
            halt(f"K={K}: r=4 alternation inequality fails")
        n_ok += 1
        worst1 = m1 if worst1 is None else min(worst1, m1)
        worst2 = m2 if worst2 is None else min(worst2, m2)
    return {"n_ok": n_ok, "worst_m1": str(worst1), "worst_m2": str(worst2)}


# ---------------------------------------------------------------------------
# Step 2c: r=5 fixed-word certificates + limit BFS + tail constants
# ---------------------------------------------------------------------------

def step2c_certs():
    w = CFG["word5"]
    rows = []
    for K in cert5_lattice():
        if K % 6 != 5:
            continue
        P = boundary_params(K)
        t_star, x0, _, _ = entry_state_exact(K)
        T = CFG["fast_orbit_T"]
        xs, it = scalar_orbit(x0, P, T)
        # itinerary tail check: last 10 B-runs form a rotation of (2,2,3,2,2)^2
        rr = [tuple(v) for v in runs_of(it)][:-1]
        tail_runs = rr[-10:]
        seq = [int(v) for ch, v in tail_runs if ch == "B"]
        if len(seq) != 5 or sorted(seq) != [2, 2, 2, 2, 3] or \
                not all(ch == "A" and v == 1 for ch, v in tail_runs
                        if ch == "A") or \
                not all(ch == "B" for ch, v in tail_runs if ch == "B"):
            halt(f"K={K}: itinerary tail runs do not match frozen r=5 word: "
                 f"{tail_runs}")
        cert = None
        n_w = len(w)
        for rot in range(n_w):
            w_rot = w[rot:] + w[:rot]
            c = certify_word(K, x0, w_rot, xs=xs, it=it, i0=0,
                             delta_frac=Fraction(*CFG["delta_frac"]))
            if c.get("ok") and c.get("t_entry") is not None:
                cert = c
                break
        if cert is None:
            halt(f"K={K}: r=5 fixed-word certificate failed (no rotation "
                 f"closed with entry)")
        if cert["t_entry"] is None or cert["t_entry"] > T - 64:
            halt(f"K={K}: r=5 entry not observed within fast window")
        rows.append({"K": K, "t_entry": cert["t_entry"],
                     "dmin": float(cert["dmin"]),
                     "lam": float(cert["lam"])})
    p = write_csv("t7_r5_fast_certs.csv", ["K", "t_entry", "dmin", "lam"],
                  rows)
    return {"rows": len(rows), "csv": p.name,
            "min_dmin": min(r["dmin"] for r in rows),
            "max_t_entry": max(r["t_entry"] for r in rows)}


# exact polynomial helpers (coefficient lists, ascending powers, Fractions)

def poly_mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def poly_add(a, b):
    n = max(len(a), len(b))
    out = [Fraction(0)] * n
    for i, x in enumerate(a):
        out[i] += x
    for i, y in enumerate(b):
        out[i] += y
    return out


def poly_scale(a, k):
    return [x * k for x in a]


def sigma5_exact_poly():
    """Generic sigma_r(n) = P(n)/Q(n) exact, r in {3,4,5}, n >= 1:
    K = 6n+r, m = 5n+r, c = (5n+r)(7n+r)/(6n+r)^2, q = (6-r)/(3(6n+r)),
    sigma = q c (1-q) / (4 (h - c q)),  h = 1/(6n+r), L = c q.
    Since h - L = h (1 - c (6-r)/3) and q = (6-r) h /3:
      sigma = c (1-q) (6-r) / (12 (1 - c (6-r)/3)).
    """
    def build(r):
        N1 = [r, 5]      # 5n+r
        N2 = [r, 7]      # 7n+r  (only correct for r <= 5: 7n+r = 2K - m)
        D = [r, 6]       # 6n+r
        num_c = poly_mul(N1, N2)
        den_c = poly_mul(D, D)
        kr = Fraction(6 - r, 3)             # (6-r)/3
        threeD = poly_scale(D, 3)
        # 1 - q = (3D - (6-r)) / (3D)
        threeD_m = list(threeD)
        threeD_m[0] -= (6 - r)
        # c(1-q) = num_c * threeD_m / (den_c * 3D)
        # 1 - c(6-r)/3 = (3 den_c - (6-r) num_c) / (3 den_c)
        # sigma = c(1-q)(6-r)/12 / (1 - c(6-r)/3)
        #        = num_c * threeD_m * (6-r) * 3 * den_c /
        #          [12 * den_c * 3D * (3 den_c - (6-r) num_c)]
        # cancel 3*den_c:
        #        = num_c * threeD_m * (6-r) /
        #          [12 * 3D * (3 den_c - (6-r) num_c)]   ... wait 3/(12)=1/4:
        # numerator: num_c * threeD_m * (6-r); denominator: 12 * D * (3den_c - (6-r)num_c)
        P = poly_scale(poly_mul(num_c, threeD_m), 6 - r)
        Q = poly_mul(poly_scale(D, 12),
                     poly_add(poly_scale(den_c, 3),
                              poly_scale(num_c, -(6 - r))))
        return P, Q
    return build(5)


SIGMA_INF = {3: Fraction(35, 4), 4: Fraction(35, 76),
             5: Fraction(35, 292)}


def sigma_poly(r):
    """sigma_r(n) = P(n)/Q(n) for class r; degree bookkeeping exact."""
    N1 = [r, 5]
    N2 = [r, 7]
    D = [r, 6]
    num_c = poly_mul(N1, N2)
    den_c = poly_mul(D, D)
    threeD = poly_scale(D, 3)
    threeD_m = list(threeD)
    threeD_m[0] -= (6 - r)
    P = poly_scale(poly_mul(num_c, threeD_m), 6 - r)
    Q = poly_mul(poly_scale(D, 12),
                 poly_add(poly_scale(den_c, 3),
                          poly_scale(num_c, -(6 - r))))
    return P, Q


def sigma5_limit():
    return SIGMA_INF[5]


def sigma5_diff_bound(n0):
    """Exact |sigma5(n) - 35/292| <= C/n for n >= n0, via polynomial coeffs."""
    P, Q = sigma_poly(5)
    inf = SIGMA_INF[5]
    den = inf.denominator
    A = poly_add(poly_scale(P, den), poly_scale(Q, -inf.numerator))
    nza = [(i, c) for i, c in enumerate(A) if c != 0]
    degA, degQ = nza[-1][0], len(Q) - 1
    assert degA == degQ - 1, (degA, degQ)
    sumA = sum(abs(c) for _, c in nza)
    C = Fraction(sumA, den * Q[-1])
    return C, A, Q


def step2d():
    """Symbolic certificates:
    (a) sigma_3 strictly decreasing and sigma_4 strictly increasing for all
        n >= 1 (polynomial sign of P(n+1)Q(n)-P(n)Q(n+1), coefficient bound
        beyond N0 + direct exact check below);
    (b) |sigma_r(n) - sigma_inf| <= C_r/n for r in {3,4,5};
    (c) consistency: sigma_poly(r) reproduces boundary_params sigma."""
    def eval_poly(A, x):
        return sum(c * x ** i for i, c in enumerate(A))

    def shift(A):
        """A(n+1) as a polynomial in n (Horner: R <- R*(n+1) + a_k)."""
        outp = [Fraction(0)]
        for k in range(len(A) - 1, -1, -1):
            outp = poly_add([Fraction(0)] + outp, outp)
            outp = poly_add(outp, [A[k]])
        return outp

    out = {}
    for r in (3, 4, 5):
        P, Q = sigma_poly(r)
        if Q[-1] < 0:
            P, Q = [-x for x in P], [-x for x in Q]
        assert Q[-1] > 0 and all(x >= 0 for x in Q), f"r={r}: Q sign"
        for n in (1, 2, 7, 40):
            K = 6 * n + r
            assert eval_poly(P, n) / eval_poly(Q, n) == \
                boundary_params(K)["sigma"], f"sigma_poly({r}) n={n}"
        # (a) monotonicity: D(n) = P(n+1)Q(n) - P(n)Q(n+1)
        if r in (3, 4):  # monotonicity needed only for these classes
            D = poly_add(poly_mul(shift(P), Q),
                         poly_scale(poly_mul(P, shift(Q)), -1))
            nz = [(i, c) for i, c in enumerate(D) if c != 0]
            lead_i, lead_c = nz[-1]
            sum_abs = sum(abs(c) for _, c in nz)
            sign_inf = 1 if lead_c > 0 else -1
            N0 = int(2 * sum_abs / abs(lead_c)) + 2
            if not all(eval_poly(D, n) * sign_inf > 0
                       for n in range(1, N0 + 1)):
                halt(f"r={r}: sigma monotonicity direct check failed "
                     f"n<=N0={N0}")
            out[f"sigma{r}_mono"] = {
                "direction": "increasing" if sign_inf > 0 else "decreasing",
                "N0": N0, "lead_coeff": str(lead_c), "deg": lead_i,
            }
        else:
            # r=5: sigma is NOT monotone; record the lattice range instead
            vals = [boundary_params(K)["sigma"] for K in range(11, 606, 6)]
            out["sigma5_range"] = {"min": str(min(vals)),
                                   "max": str(max(vals)),
                                   "limit": str(SIGMA_INF[5])}
        # (b) difference from the class limit
        inf = SIGMA_INF[r]
        den = inf.denominator
        A = poly_add(poly_scale(P, den), poly_scale(Q, -inf.numerator))
        nza = [(i, c) for i, c in enumerate(A) if c != 0]
        degA, degQ = nza[-1][0], len(Q) - 1
        sumA = sum(abs(c) for _, c in nza)
        e = degA - degQ
        assert e <= -1, f"r={r}: sigma limit mismatch (degA={degA})"
        C = Fraction(sumA, den * Q[-1])  # |sigma - inf| <= C n^e, n >= 1
        out[f"sigma{r}_diff"] = {"decay": f"n^{e}", "C": str(C),
                                 "degA": degA, "degQ": degQ}
    return out


def limit_island_bfs():
    """Exact interval covering of the limit r=5 map on the band:
    every band point reaches the certified island within CFG['bfs_depth']
    limit steps. Returns dict with closure flag and min threshold gap."""
    sig = sigma5_limit()
    beta = Fraction(1)
    return island_bfs_generic(beta, sig, "limit")


def island_bfs_generic(beta, sigma, tag, delta_frac=None, depth_cap=None):
    """Exact BFS covering of the r=5 band map with parameters (beta, sigma).
    Crossing: xi_B(m) = beta^m xi + sigma (1 + beta + ... + beta^{m-1});
    m = min{m : xi_B(m) > 1}; then A-step xi' = (3/4) xi_B(m).
    Island: certified W5 cycle intervals. delta_frac: fraction of dmin used
    as island half-width; if None, take the largest dyadic fraction in
    (0, 1/2] whose chain closes (chain = island images stay inside next
    phase interval with one-sided branch conditions)."""
    S2 = 1 + beta
    S3 = S2 + beta * beta
    band_lo = Fraction(3, 4)
    band_hi = Fraction(3, 4) * (beta + sigma)
    # m = 2 iff beta^2 xi + sigma S2 > 1; threshold:
    thr2 = (1 - sigma * S2) / (beta * beta)
    # sanity: on the band, m in {2,3}: cross by 3 for all band points,
    # and not all cross at 2 (else single cell)
    assert beta * band_lo + sigma > band_lo, "band A-step sanity"
    assert (beta ** 3) * band_lo + sigma * S3 > 1, "cross-by-3 sanity"
    assert thr2 > band_lo and thr2 < band_hi, "two-cell sanity"
    # island: cycle of word (2,2,3,2,2); compose B-steps xi -> beta xi + sigma
    seq = [2, 2, 3, 2, 2]
    lam, mu = Fraction(1), Fraction(0)
    for msteps in seq:
        for _ in range(msteps):
            lam, mu = lam * beta, mu * beta + sigma
        lam, mu = lam * Fraction(3, 4), mu * Fraction(3, 4)
    xstar = mu / (1 - lam)
    ys = [xstar]
    for msteps in seq:
        y = ys[-1]
        for _ in range(msteps):
            y = y * beta + sigma
        ys.append(Fraction(3, 4) * y)
    assert ys[5] == ys[0]
    margins = []
    for i, msteps in enumerate(seq):
        y = ys[i]
        before = y
        for j in range(msteps):
            before = before * beta + sigma
        prev = y
        for j in range(msteps - 1):
            prev = prev * beta + sigma
        margins.append(1 - prev)       # gap below 1 before last B step (>=0)
        margins.append(before - 1)     # overshoot above 1 (> 0)
        assert band_lo < y <= band_hi
    dmin = min(abs(v) for v in margins)
    assert dmin > 0

    def chain_closes(delta):
        """Does the delta-island around the W5 cycle close its 5-phase
        chain with one-sided branch conditions?"""
        for j in range(5):
            a, bb = ys[j] - delta, ys[j] + delta
            mbr = seq[j]
            good = True
            for s in range(1, mbr + 1):
                a, bb = beta * a + sigma, beta * bb + sigma
                if s < mbr:
                    if not bb <= 1:  # whole interval not yet crossed
                        good = False
                        break
                else:
                    if not a > 1:  # whole interval strictly crossed
                        good = False
                        break
            if not good:
                return False
            # A-step
            a, bb = Fraction(3, 4) * a, Fraction(3, 4) * bb
            nxt = (j + 1) % 5
            if not (ys[nxt] - delta <= a and bb <= ys[nxt] + delta):
                return False
        return True

    if delta_frac is None:
        delta_frac = Fraction(1, 2)
        while delta_frac > Fraction(1, 512):
            if chain_closes(delta_frac * dmin):
                break
            delta_frac /= 2
        else:
            halt(f"island BFS [{tag}]: no dyadic delta fraction closes")
    delta = delta_frac * dmin
    if not chain_closes(delta):
        halt(f"island BFS [{tag}]: chain does not close at {delta_frac}")
    island = [(y - delta, y + delta) for y in ys[:5]]

    def step_interval(a, b):
        """Map [a,b] through split-at-threshold B^m then A; returns images."""
        out = []
        pieces = []
        if a <= thr2 < b:
            pieces.append((a, thr2, 3))
            pieces.append((thr2, b, 2))
        elif b <= thr2:
            pieces.append((a, b, 3))
        else:
            pieces.append((a, b, 2))
        for (u, v, mbr) in pieces:
            S = sum(beta ** i for i in range(mbr))
            a2 = Fraction(3, 4) * ((beta ** mbr) * u + sigma * S)
            b2 = Fraction(3, 4) * ((beta ** mbr) * v + sigma * S)
            out.append((a2, b2))
        return out

    live = [(band_lo, band_hi, 0)]
    min_thr_gap = None
    cap = depth_cap if depth_cap is not None else CFG["bfs_depth"]
    n_pop = 0
    while live:
        lo, hi, depth = live.pop()
        n_pop += 1
        if n_pop > 200000:
            halt(f"island BFS [{tag}] population cap exceeded")
        if depth > cap:
            halt(f"island BFS [{tag}] depth exceeded at ({lo},{hi})")
        if any(l <= lo and hi <= u for l, u in island):
            continue
        gap = min(abs(lo - thr2), abs(hi - thr2))
        min_thr_gap = gap if min_thr_gap is None else min(min_thr_gap, gap)
        for a2, b2 in step_interval(lo, hi):
            live.append((a2, b2, depth + 1))
    return {"closed": True, "dmin": dmin, "delta": delta,
            "delta_frac": str(delta_frac), "min_thr_gap": min_thr_gap,
            "thr2": thr2, "xstar": xstar, "n_pop": n_pop,
            "island": [(str(l), str(u)) for l, u in island]}


def step2c_bfs():
    """Finite-K BFS covering for lattice K in [K0_5, K_tail):
    EVERY band point reaches the certified W5 island (global attraction).
    K >= K_tail is covered by the family-interval covering tail."""
    worst_gap = None
    worst_dmin = None
    n_ok = 0
    for K in cert5_lattice():
        if K_TAIL[0] is not None and K >= K_TAIL[0]:
            continue
        res = island_bfs_generic(boundary_params(K)["beta"],
                                 boundary_params(K)["sigma"], f"K={K}",
                                 depth_cap=CFG["bfs_depth_deep"])
        n_ok += 1
        worst_gap = (res["min_thr_gap"] if worst_gap is None
                     else min(worst_gap, res["min_thr_gap"]))
        worst_dmin = (res["dmin"] if worst_dmin is None
                      else min(worst_dmin, res["dmin"]))
        if n_ok % 50 == 0:
            print(f"  [2c-bfs] {n_ok} K covered (up to K={K})", flush=True)
    return {"n_ok": n_ok, "worst_min_thr_gap": str(worst_gap),
            "worst_dmin": str(worst_dmin), "K_tail": K_TAIL[0]}


def family_island_cover(K0, C5, w_frac_force=None):
    """Family-interval covering for ALL lattice K >= K0 (r = 5).

    Parameter box (proven envelope, not an assumption):
      beta_K in [1 - 1/(12 K0), 1]   (beta = 1 - q/4, q = 1/(3K))
      sigma_K in [sigma_inf - C5/n0, sigma_inf + C5/n0],  n0 = (K0-5)/6
        (|sigma_5(n) - 35/292| <= C5/n, polynomial certificate).
    The B-run threshold theta_2(beta, sigma) = (1 - sigma(1+beta))/beta^2
    is strictly decreasing in both arguments, so its family range is
    [theta_2(1, sigma_hi), theta_2(beta_lo, sigma_lo)]. Points in the
    uncertainty strip between the two extremes are evolved under BOTH
    branches (union of images); images are computed at opposite parameter
    corners (F is strictly increasing in xi, beta, sigma on the box), so
    every computed interval is a superset of every family map's image.
    Closure into the family island (limit cycle points +/- w, w the
    largest dyadic fraction of the limit dmin whose FAMILY chain closes,
    with worst-corner crossing conditions) proves: every orbit of every
    map in the box enters the island within the recorded depth and then
    follows the word W5 forever."""
    sig_inf = SIGMA_INF[5]
    n0 = (K0 - 5) // 6
    eps_sig = C5 / n0
    beta_lo = 1 - Fraction(1, 12 * K0)
    beta_hi = Fraction(1)
    sig_lo, sig_hi = sig_inf - eps_sig, sig_inf + eps_sig
    th = lambda b, s: (1 - s * (1 + b)) / (b * b)
    th_lo = th(beta_hi, sig_hi)
    th_hi = th(beta_lo, sig_lo)
    band_hi = Fraction(3, 4) * (beta_hi + sig_hi)
    # proven envelope sanity on actual lattice points
    for K in (K0, K0 + 600, K0 + 60000):
        P = boundary_params(K)
        assert K % 6 == 5
        assert beta_lo <= P["beta"] <= beta_hi
        assert sig_lo <= P["sigma"] <= sig_hi
    # limit cycle and margins
    lam, mu = Fraction(1), Fraction(0)
    for mbr in (2, 2, 3, 2, 2):
        for _ in range(mbr):
            lam, mu = lam, mu + sig_inf
        lam, mu = lam * Fraction(3, 4), mu * Fraction(3, 4)
    xstar = mu / (1 - lam)
    ys = [xstar]
    for mbr in (2, 2, 3, 2, 2):
        y = ys[-1]
        for _ in range(mbr):
            y = y + sig_inf
        ys.append(Fraction(3, 4) * y)
    assert ys[5] == ys[0]
    margins = []
    for i, mbr in enumerate((2, 2, 3, 2, 2)):
        y = ys[i]
        margins.append(1 - (y + (mbr - 1) * sig_inf))
        margins.append((y + mbr * sig_inf) - 1)
    dmin = min(abs(v) for v in margins)

    def chain_closes_family(w, cross_slack=None, incl_slack=None):
        for j in range(5):
            a, bb = ys[j] - w, ys[j] + w
            mbr = (2, 2, 3, 2, 2)[j]
            for s in range(1, mbr + 1):
                a = beta_lo * a + sig_lo        # min corner
                bb = beta_hi * bb + sig_hi      # max corner
                if s < mbr:
                    if not bb <= 1:             # no family map crosses yet
                        return False
                    if cross_slack is not None:
                        cross_slack.append(1 - bb)
                else:
                    if not a > 1:               # every family map crossed
                        return False
                    if cross_slack is not None:
                        cross_slack.append(a - 1)
            a, bb = Fraction(3, 4) * a, Fraction(3, 4) * bb
            nxt = (j + 1) % 5
            lo_n, hi_n = ys[nxt] - w, ys[nxt] + w
            if not (lo_n <= a and bb <= hi_n):
                return False
            if incl_slack is not None:
                incl_slack.append(a - lo_n)     # inclusion slack (left)
                incl_slack.append(hi_n - bb)    # inclusion slack (right)
        return True

    if w_frac_force is not None:
        frac = Fraction(*w_frac_force)
        if not chain_closes_family(frac * dmin):
            return {"ok": False, "reason": "forced w_frac chain fails"}
    else:
        frac = Fraction(1, 2)
        while frac > Fraction(1, 512):
            if chain_closes_family(frac * dmin):
                break
            frac /= 2
        else:
            return {"ok": False, "reason": "no dyadic family chain closes"}
    w = frac * dmin
    island = [(y - w, y + w) for y in ys[:5]]

    def fam_img(u, v, mbr):
        S_lo = sum(beta_lo ** i for i in range(mbr))
        S_hi = sum(beta_hi ** i for i in range(mbr))
        return (Fraction(3, 4) * (beta_lo ** mbr * u + sig_lo * S_lo),
                Fraction(3, 4) * (beta_hi ** mbr * v + sig_hi * S_hi))

    live = [(Fraction(3, 4), band_hi, 0)]
    n_pop = 0
    max_depth = 0
    while live:
        lo, hi, d = live.pop()
        n_pop += 1
        if n_pop > 200000:
            return {"ok": False, "reason": "population cap"}
        if d > CFG["bfs_depth_deep"]:
            return {"ok": False, "reason": f"depth cap at ({lo},{hi})"}
        max_depth = max(max_depth, d)
        if any(l <= lo and hi <= u for l, u in island):
            continue
        pieces = []
        if hi <= th_lo:
            pieces.append((lo, hi, 3))
        elif lo >= th_hi:
            pieces.append((lo, hi, 2))
        else:
            if lo < th_lo:
                pieces.append((lo, th_lo, 3))
            s_lo, s_hi = max(lo, th_lo), min(hi, th_hi)
            if s_lo < s_hi:  # uncertainty strip: union of both branches
                pieces.append((s_lo, s_hi, 2))
                pieces.append((s_lo, s_hi, 3))
            if hi > th_hi:
                pieces.append((th_hi, hi, 2))
        for (u, v, mbr) in pieces:
            a2, b2 = fam_img(u, v, mbr)
            live.append((a2, b2, d + 1))
    # family chain margins at the accepted half-width
    cross_slack, incl_slack = [], []
    chain_closes_family(w, cross_slack=cross_slack, incl_slack=incl_slack)
    min_cross = min(cross_slack)
    min_incl = min(incl_slack)
    return {"ok": True, "K0": K0, "w_frac": str(frac), "dmin": str(dmin),
            "th_lo": str(th_lo), "th_hi": str(th_hi),
            "strip_width": str(th_hi - th_lo), "n_pop": n_pop,
            "max_depth": max_depth,
            "min_crossing_margin": str(min_cross),
            "min_inclusion_margin": str(min_incl),
            "island": [(str(y - w), str(y + w)) for y in ys[:5]],
            "beta_lo": str(beta_lo), "sig_lo": str(sig_lo),
            "sig_hi": str(sig_hi), "band_hi": str(band_hi)}


def step2c_tail():
    """Tail layer for r=5: family-interval covering over the proven
    parameter envelope. Replaces the earlier limit-shadowing argument
    (which ignored the moving discontinuity threshold; see review P1-1).
    Exports the FIXED covering certificate (parameter box, island centers
    and half-width, thresholds, margins, depth) for independent
    verification without re-running the K0/w search."""
    C5, _, _ = sigma5_diff_bound(CFG["n0_tail"])
    K_tail = None
    res = None
    for K0 in range(65, CFG["uni_K_hi"] + 1, 6):
        if K0 % 6 != 5:
            continue
        r = family_island_cover(K0, C5)
        if r.get("ok"):
            K_tail = K0
            res = r
            break
    if K_tail is None:
        halt("family covering does not close below uni_K_hi")
    K_TAIL[0] = K_tail
    cert = {
        "K0": K_tail, "C5_per_n": str(C5), "sigma_inf": str(SIGMA_INF[5]),
        "beta_lo": res["beta_lo"], "sig_lo": res["sig_lo"],
        "sig_hi": res["sig_hi"], "th_lo": res["th_lo"],
        "th_hi": res["th_hi"], "strip_width": res["strip_width"],
        "dmin_limit": res["dmin"], "w_frac": res["w_frac"],
        "min_crossing_margin": res["min_crossing_margin"],
        "min_inclusion_margin": res["min_inclusion_margin"],
        "island": res["island"], "band_hi": res["band_hi"],
        "n_pop": res["n_pop"], "max_depth": res["max_depth"],
        "sequence": [2, 2, 3, 2, 2],
    }
    (OUT / "t7_family_cover_cert.json").write_text(
        json.dumps(cert, indent=2, ensure_ascii=False), encoding="utf-8")
    return {"C_per_n": str(C5), "K_tail": K_tail,
            "w_frac": res["w_frac"], "strip_width": res["strip_width"],
            "n_pop": res["n_pop"], "max_depth": res["max_depth"],
            "th_lo": res["th_lo"], "th_hi": res["th_hi"],
            "min_crossing_margin": res["min_crossing_margin"],
            "min_inclusion_margin": res["min_inclusion_margin"],
            "cert_json": "t7_family_cover_cert.json",
            "uni_K_hi": CFG["uni_K_hi"]}


def verify_family_cover(cert=None, hard=True):
    """Fixed-parameter re-computation check of the family-cover
    certificate.

    Scope (stated exactly): this is NOT an independent reimplementation —
    it re-runs the generator family_island_cover with the recorded (K0,
    forced w) and cross-checks EVERY proof-relevant field of the JSON.
    Crucially, the envelope constant C5 is NOT trusted from the JSON: it
    is re-derived here from the verified polynomial certificate
    (sigma5_diff_bound; whose own sign/degree checks run in step2d), and
    the recorded value must match exactly. This blocks self-consistent
    certificates built on a wrong (e.g. too small) envelope that would
    fail at large K. `cert` may be passed (dict) for negative testing;
    `hard=False` returns the check dict instead of halting."""
    if cert is None:
        cert = json.loads(
            (OUT / "t7_family_cover_cert.json").read_text(encoding="utf-8"))
    # 1. bind the envelope constant to the polynomial certificate
    C5_derived, _, _ = sigma5_diff_bound(CFG["n0_tail"])
    checks = {
        "C5_from_polynomial": cert.get("C5_per_n") == str(C5_derived),
        "sigma_inf": cert.get("sigma_inf") == str(SIGMA_INF[5]),
        "sequence": cert.get("sequence") == [2, 2, 3, 2, 2],
    }
    if not all(checks.values()):
        if hard:
            halt(f"family cover certificate verification failed: {checks}")
        return {"verified": False, "checks": checks}
    num, den = cert["w_frac"].split("/")
    r = family_island_cover(cert["K0"], C5_derived,
                            w_frac_force=(int(num), int(den)))
    # 2. cross-check every proof field against the re-computation
    checks.update({
        "ok": bool(r.get("ok")),
        "beta_lo": cert.get("beta_lo") == r.get("beta_lo"),
        "sig_lo": cert.get("sig_lo") == r.get("sig_lo"),
        "sig_hi": cert.get("sig_hi") == r.get("sig_hi"),
        "band_hi": cert.get("band_hi") == r.get("band_hi"),
        "th_lo": cert.get("th_lo") == r.get("th_lo"),
        "th_hi": cert.get("th_hi") == r.get("th_hi"),
        "strip_width": cert.get("strip_width") == r.get("strip_width"),
        "strip_consistent": (Fraction(cert["th_hi"]) -
                             Fraction(cert["th_lo"]) ==
                             Fraction(cert["strip_width"])),
        "strip_positive": Fraction(cert["th_hi"]) >
                          Fraction(cert["th_lo"]),
        "dmin_limit": cert.get("dmin_limit") == r.get("dmin"),
        "island": [list(x) for x in r.get("island", [])] ==
                  [list(x) for x in cert.get("island", [])],
        "n_pop": cert.get("n_pop") == r.get("n_pop"),
        "depth_within_cap": r.get("max_depth", 10 ** 9) <=
                            CFG["bfs_depth_deep"],
        "depth": cert.get("max_depth") == r.get("max_depth"),
        "min_crossing_margin": (
            cert.get("min_crossing_margin") ==
            r.get("min_crossing_margin") and
            Fraction(cert.get("min_crossing_margin", "-1")) > 0),
        "min_inclusion_margin": (
            cert.get("min_inclusion_margin") ==
            r.get("min_inclusion_margin") and
            Fraction(cert.get("min_inclusion_margin", "-1")) > 0),
    })
    # (the parameter box's lattice-point containment at K0, K0+600 and
    #  K0+60000 is asserted inside family_island_cover itself, hence
    #  re-checked by this re-computation; the load-bearing universal
    #  guarantee is C5_from_polynomial above)
    if not all(checks.values()):
        if hard:
            halt(f"family cover certificate verification failed: "
                 f"{ {k: v for k, v in checks.items() if not v} }")
        return {"verified": False, "checks": checks}
    return {"verified": True, "checks": checks,
            "K0": cert["K0"], "max_depth": cert["max_depth"],
            "min_inclusion_margin": cert["min_inclusion_margin"],
            "min_crossing_margin": cert["min_crossing_margin"]}


# ---------------------------------------------------------------------------
# Step 3: LOCK lattice
# ---------------------------------------------------------------------------

def step3():
    rows = []
    for K in CFG["lock_K"]:
        if K % 6 not in (0, 1, 2):
            continue
        P = boundary_params(K)
        m, h, c, q, L, a = (P["m"], P["h"], P["c"], P["q"], P["L"], P["a"])
        p = replay.exact_bin_probs(replay.near_cdf_exact, K)
        g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
        suffix0 = p[m:]
        t_lastB = None
        t_entry = None
        ok_sel_d = True
        gapprev = None
        x_mono = True
        for t in range(CFG["lock_steps"]):
            D = replay.cumsum_excess(p, K)
            s = replay.argmax_first(D)
            d = D[s - 1]
            if s in (m - 1, m) and t_entry is None:
                t_entry = t
            if s == m - 1:
                t_lastB = t
            # from max(last B, entry) + 1 onward: selection is m, d == a,
            # and |x - L| is non-increasing (x = current pre-update value)
            t_lock = t_lastB if t_lastB is not None else t_entry
            if t_lock is not None and t > t_lock:
                if s != m or d != a:
                    ok_sel_d = False
                gap = abs(p[m - 1] - L)
                if gapprev is not None and gap > gapprev:
                    x_mono = False
                gapprev = gap
            p = replay.update_state(p, g, s, ALPHA)
        # eventual fixation: last 100 selections all m
        p2 = p
        ev_m = True
        for t in range(100):
            D2 = replay.cumsum_excess(p2, K)
            if replay.argmax_first(D2) != m:
                ev_m = False
                break
            p2 = replay.update_state(p2, g, m, ALPHA)
        frozen = all(p[m + i] == suffix0[i] for i in range(K - m))
        if not (ok_sel_d and x_mono and frozen and ev_m):
            halt(f"K={K}: LOCK-side invariant failed "
                 f"(sel/d={ok_sel_d}, mono={x_mono}, frozen={frozen}, "
                 f"ev_m={ev_m})")
        rows.append({"K": K, "r": K % 6, "t_lastB": t_lastB,
                     "t_entry": t_entry, "ok": True})
    p = write_csv("t7_lock.csv", ["K", "r", "t_lastB", "t_entry", "ok"],
                  rows)
    return {"rows": len(rows), "csv": p.name}


# ---------------------------------------------------------------------------

def main():
    t0 = time.perf_counter()
    res = {"config": CFG}
    res["step0"] = step0()
    res["step1"] = step1()
    res["step2a"] = step2a()
    res["step2b"] = step2b()
    res["step2c_certs"] = step2c_certs()
    res["step2c_tail"] = step2c_tail()
    res["step2c_verify"] = verify_family_cover()
    res["step2c_bfs"] = step2c_bfs()
    res["step2d"] = step2d()
    res["step3"] = step3()
    res["events"] = EVENTS
    res["hard_halts"] = HARD
    res["elapsed_s"] = round(time.perf_counter() - t0, 1)
    res["git"] = git_sha()
    # output hashes
    res["output_sha256"] = {
        p.name: sha(p) for p in sorted(list(OUT.glob("*.csv")) +
                                       list(OUT.glob("*.json")))
        if p.name != "t7_summary.json"}
    (OUT / "t7_summary.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k != "config"},
                     indent=2, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
