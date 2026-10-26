# -*- coding: utf-8 -*-
"""v04 interval certifier v3 (review round-2 I1/I2/I3/I4/I5 rebuild).

Changes vs v2 (2b14700), mandated by ZCODE_V04_INTERVAL_REVIEW_CN:
  I1  symbolic core = symkernel (ported v0.3.0 arithmetic; squarefree
      takes the QUOTIENT; multiple/endpoint roots handled).
  I2  W0 added (U0-U2>0, U1-U2>0 on the interval); full time coverage
      recorded and PASS requires word_complete (W3 never broke early):
      {0,1}: W0, t=2: reference identity, 3..T1: W4, >=T1+1: W5
      (z <= rho^(T1+1-t0)).
  I3  explicit target/pass/open reporting; common interval scoped to
      the pass set only.
  I4  pieces: 'closed' on [lo,hi] and 'leftopen' on (lo,hi]; the tie
      family is certified as LEFT closed piece [a0-r,a0] (word = the
      smallest-index tie winner, which IS the true orbit at a0) plus
      RIGHT left-open piece (a0,a0+r]; no cross-branch continuity is
      assumed; each piece uses its own true U2.
  I5  S2 >= 0 replaced by the exact Gram identity
      (1-q)^2 (A E - d^2) == K^2 [ E2 (1-q)^2 + G2 E2 - GE^2 ],
      whose bracket is >= 0 by Cauchy-Schwarz (Gram determinant);
      an identity mismatch is a hard failure.

Records bind: engine version, K, interval + mode, word up to T1, T1,
route, fails, coverage note.  Residence: theorem A (T1'-T3').
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from symkernel import (neg_closed, neg_leftopen, nonneg_closed,  # noqa
                       nonneg_leftopen, peval, pos_closed, pos_leftopen,
                       ptrim, BudgetExceeded)
from validate_interval import PRODUCTION_ENGINE  # noqa

VERSION = PRODUCTION_ENGINE   # single source (round-3 fix 3.1)
# v3 records in interval_certs_v3.jsonl remain accepted as-is: they were
# validated by the INDEPENDENT mathematical auditor
# (tools/audit_v04_interval_v3.py, Bernstein semantics, which does NOT
# import symkernel).  v4 governs any future certificate production; the
# audited artifact SHA-256 is recorded in the review report.
ALPHA0 = F(1, 4)
KS = list(range(6, 204)) + [300, 600, 999, 1000, 1001]
RADII = [F(1, 125), F(1, 250), F(1, 500), F(1, 1000), F(1, 2000),
         F(1, 4000), F(1, 8000)]
TIE_FAMILY = {16, 32, 96, 160}       # exact s1 tie at a0 (probe range)


def Gc(x: F) -> F:
    if x <= 0:
        return F(0)
    if x < F(1, 3):
        return x / 2
    if x <= F(2, 3):
        return 3 * x / 2 - F(1, 3)
    if x <= F(5, 6):
        return 2 * x - F(2, 3)
    return F(1)


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


def step_poly(p, s, g):
    K = len(p)
    Fs = [F(0)]
    for i in range(s):
        Fs = padd(Fs, p[i])
    M = pmul([0, 1], Fs)
    out = []
    for j in range(K):
        base = p[j] if j + 1 > s else padd(p[j], pmul([0, -1], p[j]))
        out.append(padd(base, pmul(M, [g[j]])))
    return out


def U_poly(p, K, h):
    s = [F(0)]
    for i in range(K):
        dd = psub(p[i], [h])
        s = padd(s, pmul(dd, dd))
    return [K * v for v in s]


def _pos(a, lo, hi, mode):
    return pos_closed(a, lo, hi) if mode == "closed" \
        else pos_leftopen(a, lo, hi)


def _nonneg(a, lo, hi, mode):
    return nonneg_closed(a, lo, hi) if mode == "closed" \
        else nonneg_leftopen(a, lo, hi)


def _sign_pinned(a, lo, hi, mode):
    """a keeps ONE sign on the whole piece (branch pinning)."""
    if mode == "closed":
        return pos_closed(a, lo, hi) or neg_closed(a, lo, hi)
    return pos_leftopen(a, lo, hi) or neg_leftopen(a, lo, hi)


def analyze_piece(K: int, lo: F, hi: F, mode: str, rep: F, T1: int):
    """Certify ONE piece.  mode='closed': positivity on [lo,hi];
    mode='leftopen': positivity on (lo,hi].  rep = rational word
    representative inside the piece (used for all argmax decisions)."""
    h = F(1, K)
    p0 = [2 * F(1, K) - F(2 * i - 1, K * K) for i in range(1, K + 1)]
    g = [Gc(F(j, K)) - Gc(F(j - 1, K)) for j in range(1, K + 1)]
    m = -(-5 * K // 6)
    b = F(m, K)
    c = 2 * b - b * b
    q = F(6 - K % 6, 3 * K)
    gA, gB = c / (1 - q), 1 / (1 - q)
    K6 = (K == 6)
    t0 = 1 if K6 else 2
    fails = []
    # ---- W1: s0 interior; s1 pinned on the piece -------------------
    D0, acc = [], F(0)
    for j, pj in enumerate(p0):
        acc += pj
        D0.append(acc - F(j + 1, K))
    s0 = 1 + max(range(K), key=lambda j: (D0[j], -j))
    if not (1 <= s0 <= m - 2):
        return _rec(K, lo, hi, mode, T1, "OPEN", ["W1:s0-interior"],
                    None, None)
    Fs0 = sum(p0[:s0])
    ic, sl = [], []
    for j in range(K):
        Gj = sum(g[:j + 1])
        ic.append(D0[j])
        sl.append((Fs0 * Gj - F(j + 1, K) - D0[j]) if j + 1 <= s0
                  else -Fs0 * (1 - Gj))
    s1 = 1 + max(range(K), key=lambda j: (ic[j] + rep * sl[j], -j))
    if K6 and s1 not in (m - 1, m):
        return _rec(K, lo, hi, mode, T1, "OPEN", ["W1:s1-not-boundary"],
                    s1, None)
    if not K6 and not (s1 <= m - 2):
        return _rec(K, lo, hi, mode, T1, "OPEN", ["W1:s1-not-interior"],
                    s1, None)
    w1_ok = True
    for j in range(K):
        if j + 1 == s1:
            continue
        diff = [ic[s1 - 1] - ic[j], sl[s1 - 1] - sl[j]]      # affine
        if j + 1 < s1:
            ok = _pos(diff, lo, hi, mode)
        else:
            ok = _pos(diff, lo, hi, mode) or _nonneg(diff, lo, hi, mode)
        if not ok:
            w1_ok = False
            break
    if not w1_ok:
        fails.append("W1:s1-switch-inside")
    # ---- trajectory with pinned word --------------------------------
    states = [[[v] for v in p0]]
    words = [s0]
    states.append(step_poly(states[-1], s0, g))
    states.append(step_poly(states[-1], s1, g))
    words.append(s1)
    word_complete = True
    if not K6:
        p2 = states[2]
        D2, ac = [], [F(0)]
        for j in range(K):
            ac = padd(ac, p2[j])
            D2.append(psub(ac, [F(j + 1, K)]))
        s2 = 1 + max(range(K), key=lambda j: (peval(D2[j], rep), -j))
        if s2 not in (m - 1, m):
            fails.append("W2:no-entry-at-2")
        else:
            if not any(all(_pos(psub(cand, D2[j]), lo, hi, mode)
                           for j in range(m - 2))
                       for cand in (D2[m - 2], D2[m - 1])):
                fails.append("W2:entry-dominance")
    for t in range(2 if K6 else t0, T1):
        pcur = states[t]
        D, ac = [], [F(0)]
        for j in range(K):
            ac = padd(ac, pcur[j])
            D.append(psub(ac, [F(j + 1, K)]))
        st = 1 + max(range(K), key=lambda j: (peval(D[j], rep), -j))
        if st not in (m - 1, m):
            fails.append(f"W3:sel{t}-not-boundary")
            word_complete = False
            break
        diff = psub(D[m - 2], D[m - 1])
        if not _sign_pinned(diff, lo, hi, mode):
            fails.append(f"W3:branch-flip@t{t}")
            word_complete = False
            break
        states.append(step_poly(pcur, st, g))
        words.append(st)
    # ---- W0 (I2 fix) + W4 --------------------------------------------
    U2 = U_poly(states[2], K, h)
    if not _pos(psub(U_poly(states[0], K, h), U2), lo, hi, mode):
        fails.append("W0:U0<=U2")
    if not _pos(psub(U_poly(states[1], K, h), U2), lo, hi, mode):
        fails.append("W0:U1<=U2")
    if not word_complete:
        fails.append("W4:incomplete-window")
    else:
        for t in range(3, T1 + 1):
            if not _pos(psub(U_poly(states[t], K, h), U2), lo, hi, mode):
                fails.append(f"W4:U{t}<=U2")
                break
    # ---- W5 tail with I5 Gram core -----------------------------------
    pT = states[t0]
    xT = pT[m - 1]
    gam = psub([gA], pmul([gB], xT))
    e = [psub(pT[i], pmul(gam, [g[i]])) for i in range(m - 1)]
    G1 = sum(g[:m - 1])
    G2 = sum(v * v for v in g[:m - 1])
    E1 = [F(0)]
    GE = [F(0)]
    E2 = [F(0)]
    for i in range(m - 1):
        E1 = padd(E1, e[i])
        GE = padd(GE, pmul([g[i]], e[i]))
        E2 = padd(E2, pmul(e[i], e[i]))
    Sfx = K * sum((p0[i] - h) ** 2 for i in range(m, K))
    A = K * gB * gB * G2 + K
    Bc = -2 * K * gA * gB * G2 + 2 * gB * G1 - 2
    Cc = (K * gA * gA * G2 - 2 * gA * G1 + F(m - 1, K) + F(1, K) + Sfx)
    d = [K * gB * (-v) for v in GE]
    f_ = psub(pmul([K * gA], GE), E1)
    E = [K * v for v in E2]
    lhs = pmul([F((1 - q) ** 2)], psub(pmul([A], E), pmul(d, d)))
    rhs = pmul([K * K], padd(pmul([F((1 - q) ** 2)], E2),
                             psub(pmul([G2], E2), pmul(GE, GE))))
    gram_ok = (ptrim(psub(lhs, rhs)) == [])
    if not gram_ok:
        fails.append("W5:GRAM-IDENTITY-FAIL")
    AH1 = psub(pmul([2 * A], f_), pmul([Bc], d))
    S2 = psub(pmul([A], E), pmul(d, d))       # = A*H2 >= 0 by Gram
    H0 = psub(psub([Cc], U2), [F(Bc * Bc, 4 * A)])
    rho = psub([F(1)], [0, 1])
    w = [F(1)]
    for _ in range(T1 + 1 - t0):
        w = pmul(w, rho)
    w2 = pmul(w, w)
    S3 = padd(pmul([2], pmul(S2, w)), AH1)     # A * H'(w)
    Hw = padd(padd(pmul(S2, w2), pmul(AH1, w)), pmul([A], H0))
    if not _nonneg(psub([], S3), lo, hi, mode):
        fails.append("W5:S3")
    if not _pos(Hw, lo, hi, mode):
        fails.append("W5:S4")
    route = "H-convex(Gram)" if gram_ok else None
    cls = "PASS" if (not fails and word_complete) else "OPEN"
    return _rec(K, lo, hi, mode, T1, cls, sorted(set(fails)),
                words[:T1], route)


def _rec(K, lo, hi, mode, T1, cls, fails, words, route):
    return {"engine": VERSION, "K": K, "lo": str(lo), "hi": str(hi),
            "mode": mode, "T1": T1, "class": cls, "fails": fails or [],
            "word_head": words, "route": route,
            "coverage": ("t=0,1:W0; t=2:ref-identity; 3..T1:W4; "
                         ">=T1+1:W5 z<=rho^(T1+1-t0)")}


def certify_symmetric(K: int, r: F, t1_ladder):
    rec = None
    for T1 in t1_ladder:
        try:
            rec = analyze_piece(K, ALPHA0 - r, ALPHA0 + r, "closed",
                                ALPHA0, T1)
        except BudgetExceeded:
            rec = _rec(K, ALPHA0 - r, ALPHA0 + r, "closed", T1, "OPEN",
                       ["BUDGET-EXCEEDED"], None, None)
        if rec["class"] != "OPEN":
            return rec
        if not any(x.startswith(("W4", "W5")) for x in rec["fails"]):
            return rec           # structural: bigger T1 won't help
    return rec


def certify_tie_split(K: int, r: F, t1_ladder):
    """I4: LEFT closed [a0-r, a0] with left word (smallest-index tie
    winner = true orbit at a0), RIGHT left-open (a0, a0+r]."""
    left = right = None
    for T1 in t1_ladder:
        try:
            left = analyze_piece(K, ALPHA0 - r, ALPHA0, "closed",
                                 ALPHA0 - r / 2, T1)
        except BudgetExceeded:
            left = _rec(K, ALPHA0 - r, ALPHA0, "closed", T1, "OPEN",
                        ["BUDGET-EXCEEDED"], None, None)
        if left["class"] != "OPEN":
            break
        if not any(x.startswith(("W4", "W5")) for x in left["fails"]):
            break
    for T1 in t1_ladder:
        try:
            right = analyze_piece(K, ALPHA0, ALPHA0 + r, "leftopen",
                                  ALPHA0 + r / 2, T1)
        except BudgetExceeded:
            right = _rec(K, ALPHA0, ALPHA0 + r, "leftopen", T1, "OPEN",
                         ["BUDGET-EXCEEDED"], None, None)
        if right["class"] != "OPEN":
            break
        if not any(x.startswith(("W4", "W5")) for x in right["fails"]):
            break
    both = (left and right and left["class"] == "PASS"
            and right["class"] == "PASS")
    return {"engine": VERSION, "K": K, "lo": str(ALPHA0 - r),
            "hi": str(ALPHA0 + r), "mode": "tie-split",
            "T1": [left["T1"] if left else None,
                   right["T1"] if right else None],
            "class": "PASS" if both else "OPEN",
            "fails": {"left": left["fails"] if left else None,
                      "right": right["fails"] if right else None},
            "word_head": {"left": left["word_head"] if left else None,
                          "right": right["word_head"] if right else None},
            "route": "left-closed+right-leftopen",
            "coverage": ("[a0-r,a0] closed-left + (a0,a0+r] leftopen; "
                         "a0 = left word (smallest-index) + v0.3.0")}


def load_checkpoint(ckpt: Path, targets, engine=VERSION):
    """Round-5 freeze guard: structural problems CONTROL adoption.
    Returns (seen, rejected): structurally valid, on-target, current-
    engine records are adoptable as RECOVERABLE CANDIDATES; everything
    else is rejected and reported.  Mathematical acceptance remains
    the independent-audit layer."""
    seen, rejected = {}, []
    if not ckpt.exists():
        return seen, rejected
    from validate_interval import validate_interval_rows
    rows = [json.loads(line) for line in
            ckpt.read_text(encoding="utf-8").splitlines() if line.strip()]
    for rec in rows:
        probs = validate_interval_rows([rec], [rec.get("K")], engine=engine)
        if probs or rec["K"] in seen or rec["K"] not in targets:
            rejected.append((rec.get("K"),
                             probs[0] if probs else "duplicate/non-target"))
            continue
        seen[rec["K"]] = rec
    return seen, rejected


def main() -> int:
    t_start = time.time()
    # v4 production checkpoint (round-3 fix 3.3: separate file; the
    # accepted legacy v3 artifact interval_certs_v3.jsonl is preserved
    # untouched and remains validated by the INDEPENDENT auditor)
    ckpt = Path(__file__).with_name("interval_certs_v4.jsonl")
    seen, rejected = load_checkpoint(ckpt, set(KS))
    if rejected:
        print(f"resume: {len(rejected)} record(s) REJECTED and "
              f"excluded (first few): {rejected[:5]}")
    for K in KS:
        if K in seen:
            continue
        t_k = time.time()
        t1_ladder = [15, 22] if K <= 300 else [12, 15]
        rec = None
        for r in RADII:
            if K in TIE_FAMILY:
                rec = certify_tie_split(K, r, t1_ladder)
            else:
                rec = certify_symmetric(K, r, t1_ladder)
            if rec["class"] == "PASS":
                break
        with open(ckpt, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec) + "\n")
        seen[K] = rec
        print(f"K={K:>4} {rec['class']:>5} r={rec['r'] if 'r' in rec else rec['lo']} "
              f"({time.time() - t_k:.1f}s) fails={rec['fails']}", flush=True)
    # I3 coverage report (pass set only for the common interval)
    passed = sorted(K for K, rec in seen.items() if rec["class"] == "PASS")
    open_ks = sorted(set(KS) - set(passed))
    from fractions import Fraction as FF
    rmin = None
    if passed:
        rs = []
        for K in passed:
            rec = seen[K]
            rs.append(FF(rec.get("r")) if "r" in rec
                      else FF(rec["hi"]) - ALPHA0)
        rmin = min(rs)
    print(f"\ntarget set: {len(KS)} K = 6..203 all + "
          f"{[k for k in KS if k > 203]}")
    print(f"PASS set: {len(passed)} K; OPEN set ({len(open_ks)}): {open_ks}")
    if rmin is not None:
        print(f"common interval, VALID FOR PASS SET ONLY: "
              f"[1/4-{rmin}, 1/4+{rmin}] = "
              f"[{float(ALPHA0 - rmin):.7f}, {float(ALPHA0 + rmin):.7f}]")
    print(f"{time.time() - t_start:.1f}s; checkpoint: {ckpt}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
