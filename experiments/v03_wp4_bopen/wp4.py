"""

2026-10-19 status semantics (review R6): solve_word returns
status in {feasible, infeasible, error}; solver anomalies
(unbounded / verify-eq-failed / verify-le-failed) are ERROR and abort
the run — they are never counted as infeasibility. The same three-way
rule applies to the robust pass and the minimal-conflict subsystem
classification (robust_and_subsystems): subsystem masks that error are
reported separately and never count as conflict evidence;
fault_injection_check() injects each anomaly at fixed points and
asserts this routing. The fixed 1e-9
margin means infeasible = 'no orbit with all strict gaps >= 1e-9'
(strict infeasibility needs common-epsilon maximization, not done).
Counts in WP4_RESULTS.md beyond the main enumeration here (116 base
words, 32 thin-window calls) came from additional interactive
explorations of longer words — process-level evidence, reproducible
only via the saved session log; the committed entry point reproduces
the 116 + 32 core sets.
WP4 — B-open(ii): exact LP search for periodic orbits with interior selections.

Method (exact, Fractions only):
  * With a FIXED selection s the mass update is LINEAR:  p' = T_s p,
    T_s = diag(3/4 on bins<=s, 1 beyond) + (1/4) g e_s^T.  So over a word
    w = (s_0..s_{L-1}) the period map A_w = T_{s_{L-1}}..T_{s_0} is linear and
    "trajectory follows w exactly, periodically"  <=>  p0 = A_w p0 plus a
    system of AFFINE inequalities in p0 (selection order, nonnegativity,
    A1 suffix dominance, C in (b, chi)).
  * Feasibility solved exactly: nullspace of (I - A_w) (Fraction Gaussian
    elimination) -> parametrize p0 = sum xi_k v_k, one mass equality, all
    inequalities in xi; tiny two-phase Fraction simplex (Bland's rule).
  * Analytic witness: for interior j<=m-2, W_j(C) = C(1-G(y)) - (b-y) with
    y=jh; W_j linear increasing in C; W_j>0  <=>  swept updates strictly
    push D_j below the anchor a=C-b, which forbids periodic interior returns.
"""
import os, sys, random, time
from fractions import Fraction as F
from pathlib import Path

ROOT = str(Path(__file__).resolve().parents[2])
for q in (os.path.join(ROOT, "src"), os.path.join(ROOT, "experiments", "s4_w1_s1_xo_precision")):
    if q not in sys.path:
        sys.path.insert(0, q)
import replay

TOL = F(1, 10 ** 12)
ALPHA = F(1, 4)

def probs(K): return replay.exact_bin_probs(replay.xo_cdf_exact, K)
def anchor(K): return replay.anchor_m_b(K)

# ---------------- linear algebra (Fractions) ----------------
def eye(K): return [[F(1) if i == j else F(0) for j in range(K)] for i in range(K)]

def matmul(A, B):
    n, p, q = len(A), len(B), len(B[0])
    C = [[F(0)] * q for _ in range(n)]
    for i in range(n):
        Ai, Ci = A[i], C[i]
        for k in range(p):
            a = Ai[k]
            if a:
                Bk = B[k]
                for j in range(q):
                    Ci[j] += a * Bk[j]
    return C

def matvec(A, x): return [sum((a * xi for a, xi in zip(row, x)), F(0)) for row in A]

def Tsel(s, g, K):
    T = [[F(0)] * K for _ in range(K)]
    for i in range(K):
        base = F(3, 4) if i < s else F(1)
        for j in range(K):
            v = base if i == j else F(0)
            if j < s:
                v += g[i] * F(1, 4)
            T[i][j] = v
    return T

def nullspace(M):
    A = [row[:] for row in M]
    rows, cols = len(A), len(A[0])
    pivots = []; r = 0
    for c in range(cols):
        pr = None
        for i in range(r, rows):
            if A[i][c] != 0:
                pr = i; break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        pivots.append(c); r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for fc in free:
        v = [F(0)] * cols
        v[fc] = F(1)
        for i, pc in enumerate(pivots):
            v[pc] = -A[i][fc]
        basis.append(v)
    return basis

# ---------------- exact two-phase simplex (Bland) ----------------
def lp_feasible(Aeq, beq, Ale, ble, n):
    """Exact feasibility of: x>=0 (n vars), Aeq x = beq, Ale x <= ble."""
    rows = []
    for r, v in zip(Aeq, beq):
        rows.append(([F(x) for x in r], F(v), False))
    for r, v in zip(Ale, ble):
        rows.append(([F(x) for x in r], F(v), True))
    norm = []
    for r, v, isle in rows:
        if v < 0:
            r = [-x for x in r]; v = -v; sgn = F(-1)
        else:
            sgn = F(1)
        norm.append((r, v, isle, sgn))
    ns = sum(1 for _, _, isle, _ in norm if isle)
    ncols = n + ns
    tab = []; si = 0
    for r, v, isle, sgn in norm:
        row = r + [F(0)] * ns
        if isle:
            row[n + si] = sgn; si += 1
        tab.append(row + [v])
    m = len(tab)
    N = ncols + m
    # layout: [0..ncols-1] structural+slack, [ncols..ncols+m-1] artificials, [N] rhs
    T = [tab[i][:-1] + [F(1) if k == i else F(0) for k in range(m)] + [tab[i][-1]]
         for i in range(m)]
    basis = [ncols + i for i in range(m)]
    rc = [(F(1) if j >= ncols else F(0)) - sum(T[i][j] for i in range(m)) for j in range(N)]
    while True:
        ent = None
        for j in range(N):
            if rc[j] < 0:
                ent = j; break
        if ent is None:
            break
        leave = None; best = None
        for i in range(m):
            if T[i][ent] > 0:
                ratio = T[i][N] / T[i][ent]
                if best is None or ratio < best or (ratio == best and basis[i] < basis[leave]):
                    best = ratio; leave = i
        if leave is None:
            return (False, None, "unbounded")
        pv = T[leave][ent]
        T[leave] = [x / pv for x in T[leave]]
        Tl = T[leave]
        for i in range(m):
            if i != leave and T[i][ent] != 0:
                f = T[i][ent]
                T[i] = [a - f * b for a, b in zip(T[i], Tl)]
        f = rc[ent]
        rc = [a - f * b for a, b in zip(rc, Tl[:N])]
        basis[leave] = ent
    obj = F(0)
    for i in range(m):
        if basis[i] >= ncols:
            obj += T[i][N]
    if obj > 0:
        return (False, None, "infeasible")
    x = [F(0)] * n
    for i in range(m):
        if basis[i] < n:
            x[basis[i]] = T[i][N]
    for r, v in zip(Aeq, beq):
        if sum(a * xi for a, xi in zip(r, x)) != v:
            return (False, None, "verify-eq-failed")
    for r, v in zip(Ale, ble):
        if sum(a * xi for a, xi in zip(r, x)) > v:
            return (False, None, "verify-le-failed")
    return (True, x, "ok")

# ---------------- word -> LP ----------------
def word_maps(K, word):
    g = probs(K)
    Ts = {s: Tsel(s, g, K) for s in set(word)}
    P = [eye(K)]
    for s in word:
        P.append(matmul(Ts[s], P[-1]))
    A = P[len(word)]
    cum = []
    for Pt in P[:len(word)]:
        R = []; acc = [F(0)] * K
        for i in range(K):
            acc = [x + y for x, y in zip(acc, Pt[i])]
            R.append(acc[:])
        cum.append(R)
    return A, cum

def solve_word(K, word, mu=F(1, 10 ** 9), strict_all=False,
               use_sel=True, use_a1=True, use_c=True, chi=None,
               lp=lp_feasible):
    """Feasibility of a legal periodic orbit with selection word `word`.
    strict_all=False -> legal min-argmax semantics (ties allowed for j>s_t).
    Returns dict(status, p0?, why?, nullity). `lp` is injectable for the
    fault-injection self-test."""
    m, b = anchor(K); L = len(word)
    A, cum = word_maps(K, word)
    IA = [[(F(1) if i == j else F(0)) - A[i][j] for j in range(K)] for i in range(K)]
    NS = nullspace(IA)
    d = len(NS)
    if d == 0:
        return {"status": "infeasible", "why": "no-fixed-point (I-A invertible)", "nullity": 0}
    ev = [sum(v) for v in NS]
    if all(x == 0 for x in ev):
        return {"status": "infeasible", "why": "fixed points all have mass 0", "nullity": d}
    if chi is None:
        chi = F(1)
    # variables u[0..d-1], w[0..d-1] >= 0, xi = u - w
    n = 2 * d
    Aeq = [[ev[j] if j < d else -ev[j - d] for j in range(n)]]
    beq = [F(1)]
    Ale = []; ble = []
    # nonneg p
    for i in range(K):
        row = [F(0)] * n
        for k in range(d):
            row[k] = -NS[k][i]; row[d + k] = NS[k][i]
        Ale.append(row); ble.append(F(0))
    if use_sel:
        for t, s in enumerate(word):
            Rs = cum[t][s - 1]
            for j in range(1, m + 1):
                if j == s:
                    continue
                margin = mu if (strict_all or j < s) else F(0)
                Rj = cum[t][j - 1]
                dif = [a - c for a, c in zip(Rj, Rs)]
                row = [F(0)] * n
                for k in range(d):
                    c = sum(dif[i] * NS[k][i] for i in range(K))
                    row[k] = c; row[d + k] = -c
                Ale.append(row); ble.append(F(j - s, K) - margin)
    if use_a1:
        for jj in range(1, K - m + 1):
            idx = list(range(m, m + jj))
            row = [F(0)] * n
            for k in range(d):
                c = sum(NS[k][i] for i in idx)
                row[k] = c; row[d + k] = -c
            Ale.append(row); ble.append(F(jj, K) - mu)
    if use_c:
        row = [F(0)] * n
        for k in range(d):
            c = sum(NS[k][i] for i in range(m))
            row[k] = c; row[d + k] = -c
        Ale.append(row); ble.append(chi - mu)
        Ale.append([-x for x in row]); ble.append(-(b + mu))
    ok, x, why = lp(Aeq, beq, Ale, ble, n)
    if not ok:
        # status classification (2026-10-19 review): solver anomalies are
        # ERROR, not mathematical infeasibility; only a clean phase-2
        # infeasibility certificate counts as infeasible
        if why == "infeasible":
            return {"status": "infeasible", "why": why, "nullity": d}
        return {"status": "error", "why": why, "nullity": d}
    xi = [x[k] - x[d + k] for k in range(d)]
    p0 = [sum(xi[k] * NS[k][i] for k in range(d)) for i in range(K)]
    # exact verification of everything we imposed
    assert sum(p0) == 1 and all(v >= 0 for v in p0)
    Ap = matvec(A, p0)
    assert Ap == p0, "fixed point failed"
    return {"status": "feasible", "p0": p0, "nullity": d}

# ---------------- robust + subsystem classification (error-aware) ----------
def robust_and_subsystems(K, word, solver=None):
    """Robust pass + minimal-conflict subsystem classification (2026-10-19
    closeout item 2): solver statuses are three-way. `robust` is
    feasible/infeasible/error; subsystem masks with status 'infeasible' are
    conflict evidence, masks with status 'error' are recorded separately and
    are NEVER counted as conflicts. Returns (robust_tag, conflicts, errors)."""
    if solver is None:
        solver = solve_word
    rr = solver(K, list(word), mu=F(1, 10 ** 6), strict_all=True)
    if rr["status"] == "feasible":
        tag = "robust-feasible"
    elif rr["status"] == "infeasible":
        tag = "robust-infeasible"
    else:
        tag = "robust-error"   # NOT infeasibility evidence
    conflicts, errors = [], []
    for mask in range(8):
        us, ua, uc = bool(mask & 1), bool(mask & 2), bool(mask & 4)
        r2 = solver(K, list(word), mu=F(1, 10 ** 9), strict_all=False,
                    use_sel=us, use_a1=ua, use_c=uc)
        if r2["status"] == "infeasible":
            conflicts.append((us, ua, uc, r2.get("why", "ineq-system")))
        elif r2["status"] == "error":
            errors.append((us, ua, uc, r2.get("why")))
    return tag, conflicts, errors


def fault_injection_check():
    """Targeted fault-injection self-test (2026-10-19 closeout item 2).
    Injects solver anomalies at fixed points and asserts they surface as
    status='error' — never as infeasibility/conflict evidence — and that a
    silently corrupted 'feasible' answer crashes the post-solve exact
    assertions instead of passing."""
    m, b = anchor(7)
    # ground truth: pure lock word (m,) is feasible at K=7 (r=1 <= 2)
    base = solve_word(7, [m])
    assert base["status"] == "feasible", f"ground truth lost: {base}"
    # (a) lp_feasible anomalies -> solve_word must report status 'error'
    for why in ("unbounded", "verify-eq-failed", "verify-le-failed"):
        def broken(*a, _w=why, **k):
            return (False, None, _w)
        res = solve_word(7, [m], lp=broken)
        assert res["status"] == "error" and res["why"] == why, res
    # (b) corrupted 'feasible' answer -> post-solve exact assertions must fire
    def lying_lp(Aeq, beq, Ale, ble, n):
        ok, x, why = lp_feasible(Aeq, beq, Ale, ble, n)
        if ok and x is not None and len(x) > 0:
            x = list(x)
            x[0] += F(1)      # silently corrupt one variable
        return (ok, x, why)
    try:
        solve_word(7, [m], lp=lying_lp)
        raise AssertionError("fault injection: corrupted feasible answer "
                             "was NOT caught")
    except AssertionError as e:
        if "was NOT caught" in str(e):
            raise
        # any other AssertionError = the exact verification fired: expected
    # (c) robust/subsystem classification separates error from conflict
    def flaky(fail_on):
        calls = [0]
        def wrapped(K, word, **kw):
            calls[0] += 1
            if calls[0] in fail_on:
                return {"status": "error", "why": "unbounded(injected)"}
            return solve_word(K, word, **kw)
        return wrapped
    # fail the robust pass (call 1): tag must be robust-error, not infeasible
    tag, conf, errs = robust_and_subsystems(7, (m,), solver=flaky({1}))
    assert tag == "robust-error", tag
    # fail one subsystem call (calls 2..9): that mask appears in errors,
    # NOT in conflicts, and conflicts equal the clean run's conflicts
    tag0, conf0, errs0 = robust_and_subsystems(7, (m,))
    tag1, conf1, errs1 = robust_and_subsystems(7, (m,), solver=flaky({5}))
    assert tag1 == tag0, (tag0, tag1)
    assert conf1 == conf0, "error run must not gain/lose conflicts"
    assert len(errs1) == 1 and errs1[0][:3] not in [c[:3] for c in conf1]
    assert errs1[0][3] == "unbounded(injected)"
    return True


# ---------------- exact dynamics ----------------
def run_exact(p0, K, H, early_break_covered=True):
    g = probs(K); m, b = anchor(K)
    p = [F(x) for x in p0]; seq = []; stopped = None
    for t in range(H):
        D = replay.cumsum_excess(p, K)
        s = replay.argmax_first(D); d = D[s - 1]
        if d <= TOL:
            stopped = t; break
        seq.append(s)
        if early_break_covered and t >= 3 and t % 4 == 3:
            # region (R) check: every interior D_j below the conserved anchor a = D_m
            if all(D[j - 1] < D[m - 1] for j in range(1, m - 1)):
                break  # no interior selection ever again (regime where W_j>0)
        p = replay.update_state(p, g, s, ALPHA)
    return seq, stopped

def step_once(p, K, g):
    D = replay.cumsum_excess(p, K)
    s = replay.argmax_first(D)
    return replay.update_state(p, g, s, ALPHA)

def count_returns(seq, K):
    m, _ = anchor(K)
    ts = None
    for t, s in enumerate(seq):
        if s >= m - 1:
            ts = t; break
    if ts is None:
        return None, [], 0
    rets = [t for t in range(ts + 1, len(seq)) if seq[t] <= m - 2]
    ep = 0; prev = None
    for t in rets:
        if prev is None or t != prev + 1:
            ep += 1
        prev = t
    return ts, rets, ep

def fmt(v): return str(v)

# =====================================================================
if __name__ == "__main__":
    t0 = time.time()
    print("== WP4 exact LP study of B-open(ii) ==")

    # ---- self tests ----
    ok, x, why = lp_feasible([[1, 1]], [1], [], [], 2)
    assert ok and x[0] + x[1] == 1, "simplex self test 1"
    ok, x, why = lp_feasible([], [], [[1, 0]], [-1], 2)
    assert not ok, "simplex self test 2"
    ok, x, why = lp_feasible([[1, 1]], [1], [[-1, 0]], [-1], 2)
    assert ok, "simplex self test 3"
    g7 = probs(7)
    rng = random.Random(1)
    for s in (1, 3, 6, 7):
        p = [F(rng.randint(0, 50), 17) for _ in range(7)]
        p[rng.randrange(7)] += F(1)
        T = Tsel(s, g7, 7)
        assert matvec(T, p) == replay.update_state(p, g7, s, ALPHA), "linearity failed"
    fault_injection_check()
    print("[self-test] simplex + update linearity + fault injection OK")

    # ---- PART A: registered K=7 counterexample replay ----
    p0reg = [F(v, 28000) for v in (6852, 4007, 4007, 4007, 4007, 2320, 2800)]
    seq, stopped = run_exact(p0reg, 7, 60)
    ts, rets, ep = count_returns(seq, 7)
    print(f"[A] registered K=7: seq[0:12]={seq[:12]} t_*={ts} interior-return times={rets} episodes={ep}")

    # ---- PART B: analytic witness W_j ----
    print("[B] monotonicity witness W_j(C)=C(1-G(y))-(b-y), min at C->b+, j<=m-2:")
    for K in (7, 9, 11, 13, 41, 47):
        m, b = anchor(K); r = K % 6
        thr = F(6, 12 - r)
        worst = None; wj = None
        for j in range(1, m - 1):
            y = F(j, K); G = replay.xo_cdf_exact(y)
            if y <= F(2, 3):
                continue  # L_j <= 0 < a automatically
            Wb = b * (1 - G) - (b - y)
            if worst is None or Wb < worst:
                worst = Wb; wj = j
        if worst is None:
            note = "all interior y<=2/3 => L_j<=0<a automatically (witness trivial)"
        else:
            note = f"min_j W_j(b)={worst} at j={wj} (positive => interior returns die)"
        print(f"    K={K:3d} r={r} m={m} b={b} thr=6/(12-r)={thr}  b>=thr:{b>=thr}  {note}")

    # ---- PART C: LP feasibility, K in {7,9,11,13} ----
    print("[C] prescribed periodic words with interior selections, legal semantics (ties j>s allowed):")
    nfeas = 0; ninfeas = 0; nerr = 0; details = []
    for K in (7, 9, 11, 13):
        m, b = anchor(K)
        jset = sorted(set(j for j in (1, 2, m - 2) if 1 <= j <= m - 2))
        words = set()
        for j in jset:
            cand = [(j, m), (j, m - 1), (m, j), (j, j), (j, j, m), (j, m, m),
                    (j, m - 1, m), (m, j, m), (m - 1, j, m)]
            if j + 1 <= m - 2:
                cand.append((j, j + 1, m))
            for w in cand:
                if any(1 <= x <= m - 2 for x in w):
                    words.add(w)
        for w in sorted(words):
            res = solve_word(K, list(w), mu=F(1, 10 ** 9), strict_all=False)
            if res["status"] == "feasible":
                nfeas += 1
                print(f"    K={K} word={w} FEASIBLE nullity={res['nullity']} p0={res['p0']}")
                details.append((K, w, res))
            elif res["status"] == "error":
                nerr += 1
                print(f"    K={K} word={w} SOLVER ERROR [{res['why']}] "
                      f"(not counted as infeasible)")
            else:
                ninfeas += 1
        # robust pass + minimal-conflict subsystem classification on core
        # words (error-aware: robust-error and subsystem errors are reported
        # separately and abort the run — never counted as infeasibility)
        for j in jset:
            for w in [(j, m), (j, m - 1), (m, j), (j, j)]:
                tag, mins, sub_errs = robust_and_subsystems(K, w)
                if sub_errs:
                    raise SystemExit(
                        f"[C] ABORT: subsystem solver errors on K={K} "
                        f"word={w}: {sub_errs}")
                print(f"    K={K} word={w}: {tag}; minimal conflicting "
                      f"subsystems: {mins}")
    print(f"[C] total: feasible={nfeas} infeasible={ninfeas} "
          f"solver_errors={nerr}")
    if nerr:
        raise SystemExit(f"[C] ABORT: {nerr} solver errors — results "
                         f"incomplete")

    # ---- PART D: thin-window K (r=5, K>=41: b < 6/7), j = m-2 ----
    print("[D] thin-window regime (monotonicity fails for j=m-2 when C in (b, 6/7)):")
    for K in (41, 47):
        m, b = anchor(K); jstar = m - 2
        print(f"    K={K} m={m} b={b} window (b,6/7) width={F(6,7)-b}")
        for w in [(jstar, m), (m, jstar), (jstar, jstar), (jstar, m - 1), (m - 1, jstar),
                  (jstar, m, m), (m, jstar, m), (jstar, m - 1, m)]:
            for chi, chilab in ((F(6, 7), "C<6/7"), (F(1), "C<1")):
                res = solve_word(K, list(w), mu=F(1, 10 ** 9), strict_all=False, chi=chi)
                st = res["status"]
                if st == "feasible":
                    print(f"    K={K} word={w} {chilab} FEASIBLE p0={res['p0']}")
                    # verify by exact dynamics: 3+ periods
                    p0 = res["p0"]
                    seq, stopped = run_exact(p0, K, 4 * len(w) + 2, early_break_covered=False)
                    match = seq == list(w) * 4
                    print(f"      exact replay seq={seq} matches word x4: {match} stopped={stopped}")
                elif st == "error":
                    print(f"    K={K} word={w} {chilab} SOLVER ERROR "
                          f"[{res.get('why')}] — NOT infeasibility evidence")
                    raise SystemExit("[D] ABORT on solver error")
                else:
                    why = res.get("why", "ineq-system")
                    print(f"    K={K} word={w} {chilab} infeasible [{why}] nullity={res.get('nullity')}")

    # ---- PART E: targeted near-tie shelf scan ----
    # After a boundary sweep at l=m-1 (pre-sweep F_l exact), for swept interior j:
    #   D'_j - D'_l = (3/4)(D_j - D_l) + (1/4)(L_j - L_l),   L_j = F_l G_j - j h
    #   D'_j - a    = (3/4) m_j     + (1/4)(L_j - a),        m_j = D_j - a > 0
    # A return needs D_j to start ABOVE the anchor (m_j>0), stay above a
    # (m_j > theta_j := (a-L_j)/3) AND overtake l (m_j > e_l - w_j,
    # w_j := (L_j-L_l)/3). So m_j in (max(theta_j, e_l - w_j), e_l), e_l = D_l - a.
    print("[E] targeted near-tie shelf scan, K in {7,9,11,13}, exact dynamics:")

    def build_shelf(K, C, hot, el):
        """hot: list of (j, m_j) with D_j = a + m_j; other interior at a - h/3."""
        m, b = anchor(K); a = C - b; l = m - 1
        D = {0: F(0), m: a, l: a + el}
        for j, mj in hot:
            D[j] = a + mj
        for j in range(1, m):
            if j not in D:
                D[j] = a - F(1, 3 * K)
        p0 = [D[j] - D[j - 1] + F(1, K) for j in range(1, m + 1)]
        p0 += [(1 - C) / (K - m)] * (K - m)
        return p0

    best = {}
    multi = []
    nruns = 0
    import itertools as it

    def evaluate(K, C, el, hot):
        global nruns
        p0 = build_shelf(K, C, hot, el)
        if any(v < 0 for v in p0) or sum(p0) != 1:
            return None
        nruns += 1
        seq, stopped = run_exact(p0, K, 400)
        ts, rets, ep = count_returns(seq, K)
        if ts is None:
            return None
        key = (ep, len(rets))
        if K not in best or key > best[K][0]:
            best[K] = (key, p0, C, seq[:40], ts, rets)
        if ep >= 2:
            multi.append((K, C, el, hot, ep, len(rets), p0, seq[:40]))
        return key

    T = 8
    for K in (7, 9, 11, 13):
        m, b = anchor(K); l = m - 1; h = F(1, K)
        Gl = replay.xo_cdf_exact(F(l, K))
        hits = {}
        for C in [b + (1 - b) * F(k, 8) for k in range(1, 8)]:
            a = C - b
            for el in [h * F(u, 16) for u in (2, 4, 6, 8, 10, 12, 14)]:
                Fl = C - h + el  # p_m = h - el
                Ll = Fl * Gl - l * h
                info = {}
                for j in range(1, m - 1):
                    Gj = replay.xo_cdf_exact(F(j, K))
                    Lj = Fl * Gj - j * h
                    th = (a - Lj) / 3
                    # return possible iff D'_j > a  <=>  m_j > theta_j, with m_j < e_l
                    if 0 <= th < el:
                        info[j] = th
                if not info:
                    continue
                cand = sorted(info, key=lambda j: el - info[j])[-3:]
                cfg_hit = []
                for j in cand:
                    lo = info[j]
                    for u in range(1, 7):
                        mj = lo + (el - lo) * F(u, T)
                        key = evaluate(K, C, el, [(j, mj)])
                        if key and key[0] >= 1:
                            cfg_hit.append((j, mj))
                if cfg_hit:
                    hits[(C, el)] = (info, cfg_hit)
        # phase 2: pairs/triples among coordinates that actually returned
        for (C, el), (info, cfg_hit) in hits.items():
            js = sorted(set(j for j, _ in cfg_hit))
            ms = sorted(set(mj for _, mj in cfg_hit))
            for j1, j2 in it.combinations(js, 2):
                for m1 in ms:
                    for u2 in range(1, 7):
                        lo2 = info[j2]
                        m2 = lo2 + (el - lo2) * F(u2, T)
                        evaluate(K, C, el, [(j1, m1), (j2, m2)])
                        evaluate(K, C, el, [(j2, m1), (j1, m2)])
            if len(js) >= 3:
                for combo in it.combinations(js, 3):
                    for us in it.product((2, 5), repeat=3):
                        hot = [(j, info[j] + (el - info[j]) * F(u, T))
                               for j, u in zip(combo, us)]
                        evaluate(K, C, el, hot)
        if K in best:
            (ep, nret), p0, C, seqh, ts, rets = best[K]
            print(f"    K={K}: best episodes={ep} total interior returns={nret} at C={C}"
                  f" t_*={ts} return times={rets[:12]}")
            print(f"      p0={p0}")
            print(f"      seq[:30]={seqh[:30]}")
        else:
            print(f"    K={K}: no boundary entry observed")
    print(f"    [scan] runs={nruns} states with >=2 return-episodes: {len(multi)}")
    for (K, C, el, hot, ep, nret, p0, seqh) in multi[:6]:
        print(f"    MULTI K={K} C={C} el={el} hot={hot} episodes={ep} total={nret}")
        print(f"      p0={p0}")
        print(f"      seq[:30]={seqh[:30]}")
    # champion long runs
    print("[E2] champions re-run 2000 steps:")
    for K, v in sorted(best.items()):
        (ep, nret), p0, C, _, _, _ = v
        seq, stopped = run_exact(p0, K, 2000)
        ts, rets, ep2 = count_returns(seq, K)
        print(f"    K={K} C={C}: 2000-step: t_*={ts} episodes={ep2} total returns={len(rets)}"
              f" stopped={stopped} return times={rets[:15]}")
    print(f"[wall] {time.time()-t0:.1f}s")
