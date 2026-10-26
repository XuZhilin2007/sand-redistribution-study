# -*- coding: utf-8 -*-
"""v04-B: general-alpha strengthened D-decomp — exact branch-by-branch
verification, direct back-substitution of the quadratic identity, and the
generalised square-completion stopping certificate prototype.

CLAIM UNDER TEST (derived by hand, verified here against true orbits):
  Canonical XO, K>=6, alpha in (0,1), exact arithmetic, smallest-index
  argmax. Invariants (I) [s_t<=m, F_{t,m}=c, D_{t,m}=a] and region (R)
  are alpha-free; hence entry into {l,m}-only selection occurs in finite
  time t*. From any t0 with two-branch residence:
      p_{t,i} = gamma(x_t) g_i + (1-alpha)^{t-t0} e_i   (i < m),
      gamma(x) = (c-x)/(1-q),  e_i = p_{t0,i} - gamma(x_{t0}) g_i,
  suffix i>m frozen at p0, and
      U(p_t) = Q(x_t) + 2 z Phi(x_t) + K E2 z^2,  z=(1-alpha)^{t-t0},
  with alpha-FREE coefficient assembly (Q quadratic, Phi affine).

Stopping certificate prototype (generalised S1-S4, anchor t0):
  S1 A>0; S2 H2=E-d^2/A>0; S3 H'(1-alpha)<=0; S4 H(1-alpha)>=0;
  plus prefix dominance U(t)>=U(t0) for t<t0.
  => U(p_t)>=U(p_{t0}) for ALL t>=t0 (stop-at-t0 globally optimal).

This script verifies each layer exactly on true orbits and reports where
the certificate holds / fails.  Exit 0 = all layers verified as claimed
(failures of S-conditions are REPORTED, not errors).
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path

T_WINDOW = 120
CASES = [(K, a) for K in (9, 10, 11, 12, 29, 65)
         for a in (F(1, 8), F(1, 5), F(1, 4), F(1, 3), F(2, 5), F(1, 2))]


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


def run(K: int, alpha: F, T: int = T_WINDOW):
    p0 = [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]
    g = [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]
    m = -(-5 * K // 6)
    h = F(1, K)
    b = F(m, K)
    c = 2 * b - b * b
    q = F(6 - K % 6, 3 * K)
    states = [p0]
    sels = []
    for _ in range(T):
        p = states[-1]
        D, acc, best = [], F(0), 0
        for j, pj in enumerate(p):
            acc += pj
            d = acc - F(j + 1, K)
            D.append(d)
            if d > D[best]:
                best = j
        s = best + 1
        F_s = sum(p[:s])
        M = alpha * F_s
        states.append([((1 - alpha) * pj if j < s else pj) + M * gj
                       for j, (pj, gj) in enumerate(zip(p, g))])
        sels.append(s)
    return states, sels, m, h, c, q


def analyze(K: int, alpha: F):
    states, sels, m, h, c, q = run(K, alpha)
    T = len(states) - 1
    U = [K * sum((y - h) ** 2 for y in st) for st in states]
    # entry + residence check
    t_star = next((t for t, s in enumerate(sels) if s in (m - 1, m)), None)
    residence = t_star is not None and all(
        s in (m - 1, m) for s in sels[t_star:])
    out = {"K": K, "alpha": f"{alpha.numerator}/{alpha.denominator}",
           "t_star": t_star, "residence": residence, "layers": []}
    if t_star is None:
        out["layers"].append("NO-ENTRY-IN-WINDOW")
        return out
    rho = 1 - alpha
    gam = lambda x: (c - x) / (1 - q)
    p0 = states[0]
    g = None  # recompute g
    gg = [G(F(j, K)) - G(F(j - 1, K)) for j in range(1, K + 1)]
    # layer 1: decomp componentwise from t0 = t_star
    e = [states[t_star][i] - gam(states[t_star][m - 1]) * gg[i]
         for i in range(m - 1)]
    dec_ok = True
    for t in range(t_star, min(t_star + 25, T + 1)):
        x = states[t][m - 1]
        z = rho ** (t - t_star)
        for i in range(m - 1):
            if states[t][i] != gam(x) * gg[i] + z * e[i]:
                dec_ok = False
                out["layers"].append(
                    f"DECOMP-FAIL t={t} i={i}")
                break
        if not dec_ok:
            break
    # layer 2: identity back-substitution
    G1 = sum(gg[:m - 1])
    G2 = sum(v * v for v in gg[:m - 1])
    E1 = sum(e)
    GE = sum(gg[i] * e[i] for i in range(m - 1))
    E2 = sum(v * v for v in e)
    Sfx = K * sum((p0[i] - h) ** 2 for i in range(m, K))
    gA, gB = c / (1 - q), 1 / (1 - q)
    A = K * gB * gB * G2 + K
    B = -2 * K * gA * gB * G2 + 2 * gB * G1 - 2
    Cc = (K * gA * gA * G2 - 2 * gA * G1 + F(m - 1, K) + F(1, K) + Sfx)
    d = -K * gB * GE
    f_ = K * gA * GE - E1
    E = K * E2
    id_ok = True
    for t in range(t_star, min(t_star + 8, T + 1)):
        x = states[t][m - 1]
        z = rho ** (t - t_star)
        quad = (A * x * x + B * x + Cc) + 2 * z * (d * x + f_) + E * z * z
        if quad != U[t]:
            id_ok = False
            out["layers"].append(f"IDENTITY-FAIL t={t}")
            break
    # layer 3: generalised square-completion certificate, anchor t0=t_star
    H2 = E - d * d / A
    H1 = 2 * f_ - B * d / A
    H0 = Cc - U[t_star] - B * B / (4 * A)
    Hp = lambda z: 2 * H2 * z + H1
    S = {"S1_A>0": A > 0, "S2_H2>0": H2 > 0,
         "S3_H'(1-a)<=0": Hp(rho) <= 0, "S4_H(1-a)>=0":
         H2 * rho * rho + H1 * rho + H0 >= 0}
    prefix_ok = all(U[t] >= U[t_star] for t in range(t_star))
    # empirical: does U dip below U(t0) after t0 within window?
    dip = min(range(t_star, T + 1), key=lambda t: U[t])
    out.update({
        "decomp_ok": dec_ok, "identity_ok": id_ok,
        "cert_S": {k: bool(v) for k, v in S.items()},
        "cert_all": all(S.values()) and prefix_ok,
        "prefix_ok": bool(prefix_ok),
        "post_t0_dip_t": dip,
        "post_t0_dip_below": U[dip] < U[t_star],
    })
    return out


def count_failures(results):
    """Acceptance gate over structured boolean fields ONLY (R1 repair).

    Implementation failures: decomp_ok / identity_ok False (message text
    is never parsed).  Logical contradiction: certificate claims
    U(t)>=U(t0) for all t (cert_all) yet the window shows a real dip
    below U(t0) -- hard failure.  S-conditions merely failing open
    (cert_all False without implementation failure) is a legitimate
    outcome, NOT counted here.
    """
    hard = []
    for r in results:
        if r.get("t_star") is None:
            continue        # no entry in window: limitation, not failure
        if not r.get("decomp_ok", False):
            hard.append((r["K"], r["alpha"], "DECOMP"))
        if not r.get("identity_ok", False):
            hard.append((r["K"], r["alpha"], "IDENTITY"))
        if r.get("cert_all") and r.get("post_t0_dip_below"):
            hard.append((r["K"], r["alpha"], "CERT-CONTRADICTION"))
        if r.get("layers") and any(
                l.startswith(("DECOMP-FAIL", "IDENTITY-FAIL"))
                for l in r["layers"]):
            hard.append((r["K"], r["alpha"], "LAYER-FLAGGED"))
    return hard


def main() -> int:
    results = [analyze(K, a) for K, a in CASES]
    for r in results:
        print(f"K={r['K']:>3} a={r['alpha']:>4} t*={r['t_star']} "
              f"res={r['residence']} dec={r.get('decomp_ok')} "
              f"id={r.get('identity_ok')} "
              f"cert={r.get('cert_all')} dipT={r.get('post_t0_dip_t')} "
              f"dip<{r.get('post_t0_dip_below')}")
    n_cert = sum(1 for r in results if r.get("cert_all"))
    dest = Path(__file__).with_name("decomp_alpha.json")
    dest.write_text(json.dumps(results, indent=1), encoding="utf-8")
    print("wrote", dest)
    print(f"cases={len(results)} cert_all=True: {n_cert}")
    hard = count_failures(results)
    print("HARD FAILURES (decomp/identity/contradiction):", len(hard), hard)
    return 1 if hard else 0


if __name__ == "__main__":
    sys.exit(main())
