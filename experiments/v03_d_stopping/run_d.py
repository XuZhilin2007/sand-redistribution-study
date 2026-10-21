"""v0.3 Theorem D runner (work item 4): exact stopping-rule analysis.

STATUS (2026-10-20, review 9): this runner holds the per-lattice-point
certificates (203 K, all-time coverage).  The LOCK half of the all-K
statement is CLOSED by Theorem B-lock (docs/V03_B_LOCK_ALLK.md, real
orbit, 96 sub-family certificate sets); the CYCLE half remains a
finite-lattice certification here -- all-K is open.  "t_off = 2" means
U(p_t) >= U(p_2) for all t; it does NOT assert U monotone in time on
CYCLE orbits.  LOCK closed-form decay uses the start of PERMANENT
m-selection (B-lock's t_lock), which differs from the first boundary
entry t* (e.g. K = 6: t* = 1 but t_off = 2).

Uniformity measure U(p) = K * sum_i (p_i - 1/K)^2 (the repo's U_density).

DISCOVERY being verified exactly on the lattice (candidate result D):
  for every canonical K, the uniformity-optimal stopping time is t = 2
  (which equals the boundary-entry time t* for K >= 7), and the
  asymptotic regime is 2.5-2.9x worse.

Checks per K (exact Fractions):
  D1 (LOCK): U(t) decomposition identity — with tau = t - t* - 1,
      p_t = p_inf + rho^tau * e (rho = 3/4, e = p_{t*+1} - p_inf),
      U(t) - U_inf = 2 rho^tau S' + rho^{2 tau} E',
      S' = K sum e_i (p_inf_i - h), E' = K sum e_i^2. Verify the identity
      on several taus; record sign(S'); certify
      U(2) <= U_inf - S'^2/E' (quadratic branch certificate) and
      U(2) = min_{0<=t<=300} U(t).
  D2 (CYCLE): exact argmin over [0, T]; asymptotic cycle minimum at ALL
      n phases via the per-phase two-epoch linear solve
      (D_lim = D2 + (D2 - D1)/(1 - (3/4)^n)); certify U(2) <= cycle-min
      and U(2) = min over the window.
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import time
from fractions import Fraction as F
import sys as _sys
try:
    _sys.set_int_max_str_digits(1_000_000)
except (AttributeError, ValueError):
    pass
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "s4_w1_s1_xo_precision"))
sys.path.insert(0, str(REPO / "experiments" / "v03_t7_cyclicity"))

import replay  # noqa: E402
from sand_m0.word_cyclicity import (  # noqa: E402
    ALPHA, boundary_params, compose_word, detect_eventual_word,
    scalar_orbit,
)

CENSUS_CSV = (REPO / "experiments" / "v03_t7_cyclicity" / "results" /
              "t7_census_certificates.csv")


def _census_word(K):
    """Committed T7 census certificate word for K (None if absent)."""
    import csv
    if not CENSUS_CSV.exists():
        return None
    for row in csv.DictReader(CENSUS_CSV.open(encoding="utf-8")):
        if int(row["K"]) == K:
            return row["word"]
    return None


def _same_rotation(w1, w2):
    if len(w1) != len(w2):
        return False
    return any(w1[i:] + w1[:i] == w2 for i in range(len(w1)))

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
CFG = {
    "K_lo": 6, "K_hi": 203, "anchors": [300, 600, 999, 1000, 1001],
    "T_lock": 300, "T_cycle": 400, "warmup_cycles": 40,
}
RHO = F(3, 4)


def U(p, K):
    h = F(1, K)
    return K * sum((pi - h) ** 2 for pi in p)


def trajectory(K, T):
    P = boundary_params(K)
    m = P["m"]
    p = replay.exact_bin_probs(replay.near_cdf_exact, K)
    g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
    states = [list(p)]
    t_star = None
    for t in range(T):
        D = replay.cumsum_excess(p, K)
        s = replay.argmax_first(D)
        if t_star is None and s in (m - 1, m):
            t_star = t
        if D[s - 1] <= F(1, 10 ** 12):
            break
        p = replay.update_state(p, g, s, ALPHA)
        states.append(list(p))
    return states, t_star, P, g


def analyze_lock(K):
    P = boundary_params(K)
    m, h, c = P["m"], P["h"], P["c"]
    states, t_star, P2, g = trajectory(K, CFG["T_lock"])
    assert t_star is not None
    p_inf = [c * gi for gi in g[:m]] + list(states[0][m:])
    U_inf = U(p_inf, K)
    # e and the identity
    p1 = states[t_star + 1]
    e = [p1[i] - p_inf[i] for i in range(K)]
    Sp = K * sum(e[i] * (p_inf[i] - h) for i in range(K))
    Ep = K * sum(ei * ei for ei in e)
    # verify identity at several taus
    for tau in (1, 5, 25, 100):
        t = t_star + 1 + tau
        if t >= len(states):
            break
        lhs = U(states[t], K) - U_inf
        rhs = 2 * RHO ** tau * Sp + RHO ** (2 * tau) * Ep
        assert lhs == rhs, f"identity fails K={K} tau={tau}"
    Uvals = [U(st, K) for st in states]
    best_t = min(range(len(Uvals)), key=lambda t: Uvals[t])
    # exact tail certificate: for tau beyond the window, U - U_inf equals
    # q(u) = 2 S' u + E' u^2 with u = rho^tau <= rho^Twindow; the minimum
    # of q over (0, u_cap] is attained at u_cap or at the vertex
    # u* = -S'/E' if u* <= u_cap. Certify U(2) <= U_inf + that minimum.
    u_cap = RHO ** (len(states) - 1 - t_star - 1)
    cands = [2 * Sp * u_cap + Ep * u_cap * u_cap]
    u_star = -Sp / Ep if Ep else None
    if u_star is not None and 0 < u_star <= u_cap:
        cands.append(-Sp * Sp / Ep)
    min_tail = min(cands)
    cert = (Uvals[2] <= U_inf + min_tail and best_t == 2)
    return {"K": K, "class": "LOCK", "t_star": t_star, "best_t": best_t,
            "U2_over_Uinf": float(Uvals[2] / U_inf), "S_sign":
            "neg" if Sp < 0 else ("pos" if Sp > 0 else "zero"),
            "quad_cert": cert, "n_states": len(states)}


def analyze_cycle(K):
    """Exact CYCLE analysis.

    The asymptotic periodic profile is computed EXACTLY per word phase by
    the geometric-sum solve for the always-swept interior coordinates
    (each contracts by exactly 3/4 per step under the exactly-periodic
    driving reconstructed from the certified periodic boundary orbit),
    then BACK-CHECKED by applying the word's exact update to the rebuilt
    profile (must reproduce the next phase exactly). The FOREVER-license
    that the actual trajectory follows the word with a fixed phase from
    t0 on is a RE-DERIVED interval-invariance certificate (T7-style,
    self-contained here): the word's delta-island is rebuilt and the
    trajectory's boundary coordinate at t0 must lie inside it at the
    aligned phase. The infinite tail is covered by a rigorous
    triangular-decay envelope (see the derivation comment at the
    envelope block): interior differences decay by (3/4)^dt, the
    boundary difference by mu^dt with mu = max(3/4, beta) < 1, and the
    one-way coupling (boundary -> interior, |coefficient| <= 1/4 per
    B-step) contributes (dx0/(4 mu)) * dt * mu^dt, so the envelope is
    g(dt) = (alpha + beta_c * dt) * mu^dt with alpha = dint0 + dx0 and
    beta_c = dx0/(4 mu), whose unimodal global maximum is found by exact
    scanning."""
    P = boundary_params(K)
    m, h, c, q = P["m"], P["h"], P["c"], P["q"]
    states, t_star, P2, g = trajectory(K, CFG["T_cycle"])
    assert t_star is not None
    x0 = states[t_star][m - 1]
    xs, it = scalar_orbit(x0, P, 2400)
    i0, n, w = detect_eventual_word(it, rep=12)
    assert w is not None, f"no word K={K}"
    # license: itinerary periodicity forever is the committed T7 census
    # certificate; verify the committed word matches (rotation-insensitive)
    census = _census_word(K)
    if census is not None:
        assert _same_rotation(w, census), \
            f"K={K}: word differs from committed T7 census certificate"
    # exact periodic boundary orbit (word affine fixed point + orbit)
    lam_w, mu_w = compose_word(w, P)
    xstar = mu_w / (1 - lam_w)
    ps = [xstar]
    for br in w:
        if br == "A":
            ps.append(F(3, 4) * ps[-1] + P["L"] * F(1, 4))
        else:
            ps.append(ps[-1] + (c - ps[-1]) * q * F(1, 4))
    assert ps[n] == ps[0]
    # periodic driving per word position: F_s = c (A) or c - ps[j] (B)
    G = [sum(g[:j + 1]) for j in range(K)]

    def b_of(jpos, col):
        s_m = (w[jpos] == "A")
        Fs = c if s_m else c - ps[jpos]
        return F(1, 4) * (Fs * G[col] - (col + 1) * h)

    rn = RHO ** n
    # exact periodic profile per word position
    profiles = []
    for jpos in range(n):
        Dper = [None] * (m - 1)  # interior coordinates 0..m-2 (cols j<=m-2)
        for col in range(m - 1):
            tot = F(0)
            for i in range(n):
                tot += RHO ** i * b_of((jpos - 1 - i) % n, col)
            Dper[col] = tot / (1 - rn)
        # suffix j > m-1 uses frozen canonical values; reconstruct p
        p_prof = [None] * K
        prevF = F(0)
        for col in range(m - 1):
            Fj = Dper[col] + (col + 1) * h
            p_prof[col] = Fj - prevF
            prevF = Fj
        p_prof[m - 1] = ps[jpos]           # boundary bin
        prevF = prevF + ps[jpos]
        for col in range(m, K):            # frozen suffix
            p_prof[col] = states[0][col]
        profiles.append(p_prof)
    # BACK-CHECK: word step on profile jpos -> profile jpos+1, exactly
    for jpos in range(n):
        s_sel = m if w[jpos] == "A" else m - 1
        nxt = replay.update_state(profiles[jpos], g, s_sel, ALPHA)
        tgt = profiles[(jpos + 1) % n]
        assert nxt == tgt, f"periodic-profile back-check fails K={K} j={jpos}"
    U_cycle = [U(prof, K) for prof in profiles]
    cyc_min = min(U_cycle)
    # trajectory alignment: find word offset at the window end t0
    Uvals = [U(st, K) for st in states]
    best_t = min(range(len(Uvals)), key=lambda t: Uvals[t])
    t0 = len(states) - 1 - 2 * n
    sel_tail = [replay.argmax_first(replay.cumsum_excess(states[t], K))
                for t in range(t0, len(states))]
    offs = None
    for cand in range(n):
        if all(sel_tail[i] == (m if w[(cand + i) % n] == "A" else m - 1)
               for i in range(len(sel_tail))):
            offs = cand
            break
    assert offs is not None, f"K={K}: trajectory tail does not follow w"
    # FOREVER-LICENSE (self-contained T7-style interval certificate):
    # rebuild the word's delta-island; if the trajectory's boundary
    # coordinate at t0 lies inside it at the aligned phase, the word and
    # phase hold for all t >= t0 (island chain invariance).
    dmin_w = min(abs(pp - h) for pp in ps[:n])
    delta_w = dmin_w / 64
    x_t0 = states[t0][m - 1]
    license_ok = abs(x_t0 - ps[offs]) <= delta_w
    assert license_ok, \
        f"K={K}: t0 boundary coordinate outside word island (license fails)"
    # exact deviations at t0
    q0 = profiles[offs]
    dx0 = abs(x_t0 - q0[m - 1])
    dint0 = max((abs(states[t0][col] - q0[col]) for col in range(m - 1)),
                default=F(0))
    # Derivation of the rigorous triangular tail envelope (dt = t - t0).
    # Both the trajectory and the periodic profile follow the SAME word
    # from t0 on (license above), so phase-aligned deviations satisfy:
    #   |x-dev|(dt)   <= dx0 * mu^dt,
    #       mu = max(3/4, beta) < 1   (every branch multiplies the
    #       boundary deviation by 3/4 or beta)
    #   |int-dev|(dt) <= (3/4)^dt * dint0
    #       + sum_{s<dt} (3/4)^{dt-1-s} * (1/4) * dx0 * mu^s
    #     (interior coordinates contract by exactly 3/4 per step; a B-step
    #     injects |d D'_j / dx| = G_j/4 <= 1/4 of the boundary deviation;
    #     every summand satisfies (3/4)^{dt-1-s} mu^s <= mu^{dt-1} since
    #     3/4 <= mu, hence the sum is <= (dx0/4) * dt * mu^{dt-1})
    #   so  |inf-norm deviation|(dt) <= g(dt) with the CONSERVATIVE
    #   coefficient (dx0/(4 mu)) absorbing the mu^{-1}:
    #       g(dt) = (alpha + beta_c * dt) * mu^dt,
    #       alpha = dint0 + dx0,  beta_c = dx0 / (4 mu).
    # g is unimodal (log-concave factor times linear); its global maximum
    # over dt >= 0 is found by exact increase-while-growing scanning.
    mu = max(RHO, P["beta"])
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
    sq_ceil = math.isqrt(K - 1) + 1        # ceil(sqrt(K)), integer
    D2 = F(sq_ceil) * g_max                # ||p_t - q_phase||_2, all t>=t0
    # Uniform lower bound on U via inner-product expansion and
    # Cauchy-Schwarz (no reverse-triangle squaring):
    #   U(p_t) = K*||q - h + d||^2
    #          = K*(Q2 + 2<q-h, d> + ||d||^2)
    #         >= K*(Q2 - 2*||q-h||*||d||)          (Cauchy-Schwarz)
    #         >= K*(Q2 - (1 + Q2)*D2)              (2*sqrt(Q2) <= 1 + Q2)
    # with Q2 = ||q - h||^2 per phase and ||d|| <= D2.
    tail_lb = None
    for jpos in range(n):
        Q2 = sum((qi - h) ** 2 for qi in profiles[jpos])
        lb = K * Q2 - K * (1 + Q2) * D2
        tail_lb = lb if tail_lb is None else min(tail_lb, lb)
    cert = (Uvals[2] <= tail_lb and best_t == 2)
    # ratio reported against the exact cycle minimum
    return {"K": K, "class": "CYCLE", "t_star": t_star, "n": n,
            "word": w, "xstar": str(xstar),
            "phase_offset": offs, "t0": t0,
            "island_delta": str(delta_w),
            "license_margin": str(delta_w - abs(x_t0 - ps[offs])),
            "best_t": best_t, "U2_over_cycmin": float(Uvals[2] / cyc_min),
            "lam_word": str(lam_w), "mu_env": str(mu),
            "g_max": float(g_max), "dint0": float(dint0), "dx0": float(dx0),
            "tail_lb": str(tail_lb), "U2": str(Uvals[2]),
            "cycle_cert": cert, "n_states": len(Uvals)}


def main():
    t0 = time.perf_counter()
    rows = []
    cert_json = []
    for K in list(range(CFG["K_lo"], CFG["K_hi"] + 1)) + CFG["anchors"]:
        r = analyze_lock(K) if K % 6 in (0, 1, 2) else analyze_cycle(K)
        rows.append(r)
        # certificate record: exact evidence for independent verification
        rec = {"K": K, "class": r["class"], "best_t": r["best_t"],
               "t_star": r["t_star"]}
        if r["class"] == "CYCLE":
            rec.update({k: r[k] for k in
                        ("n", "word", "xstar", "phase_offset", "t0",
                         "island_delta", "license_margin", "lam_word",
                         "mu_env", "dint0", "dx0", "tail_lb", "U2")})
        else:
            rec["t_lock_quad"] = True
        cert_json.append(rec)
        ok = r.get("quad_cert") or r.get("cycle_cert")
        if not ok or r["best_t"] != 2:
            print(f"HARD FAIL K={K}: {r}", flush=True)
            return 1
    fields = sorted({k for r in rows for k in r})
    with (OUT / "d_rows.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    (OUT / "d_certificates.json").write_text(
        json.dumps(cert_json, indent=1, ensure_ascii=False), encoding="utf-8")
    res = {
        "rows": len(rows), "all_best_t_2": all(r["best_t"] == 2 for r in rows),
        "all_certified": all(r.get("quad_cert") or r.get("cycle_cert")
                             for r in rows),
        "U2_ratio_min": min(float(r.get("U2_over_Uinf",
                                        r.get("U2_over_cycmin")))
                            for r in rows),
        "U2_ratio_max": max(float(r.get("U2_over_Uinf",
                                        r.get("U2_over_cycmin")))
                            for r in rows),
        "elapsed_s": round(time.perf_counter() - t0, 1),
        "git": subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                              capture_output=True, text=True).stdout.strip(),
        "source_files_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(__file__)]},
        "outputs": ["d_rows.csv", "d_certificates.json", "d_summary.json"],
    }
    (OUT / "d_summary.json").write_text(
        json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res))
    return 0 if res["all_best_t_2"] and res["all_certified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
