"""v0.3 frontier tests: word-cyclicity module, T7 certificates, Theorem B
initial-state machinery, Theorem C bounds/classifiers. Exact arithmetic only.

Run:  python tests/test_v03_frontier.py   (prints [PASS]/[FAIL] lines)
"""
from __future__ import annotations

import csv
import sys
from fractions import Fraction as F
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "s4_w1_s1_xo_precision"))
sys.path.insert(0, str(REPO / "experiments" / "v03_t7_cyclicity"))
sys.path.insert(0, str(REPO / "experiments" / "v03_b_robustness"))
sys.path.insert(0, str(REPO / "experiments" / "v03_c_alpha_window"))
sys.path.insert(0, str(REPO / "experiments" / "v03_d_stopping"))

import replay  # noqa: E402
from sand_m0.word_cyclicity import (  # noqa: E402
    ALPHA, boundary_params, certify_word, compose_word,
    detect_eventual_word, pure_window_membership, runs_of, scalar_orbit,
    scalar_step, sigma_strictly_monotone,
)


def check(name, cond):
    print(("[PASS]" if cond else "[FAIL]") + f" {name}")
    return cond


def t_params():
    P = boundary_params(9)
    ok = True
    ok &= check("boundary_params K=9 m", P["m"] == 8)
    ok &= check("boundary_params K=9 c=80/81", P["c"] == F(80, 81))
    ok &= check("boundary_params K=9 q=1/9", P["q"] == F(1, 9))
    ok &= check("boundary_params K=9 sigma=160/9", P["sigma"] == F(160, 9))
    ok &= check("boundary_params K=10 sigma=77/170",
                boundary_params(10)["sigma"] == F(77, 170))
    return ok


def t_scalar():
    P = boundary_params(9)
    x1, b1 = scalar_step(F(1, 10), P)   # x < h -> B
    x2, b2 = scalar_step(F(1, 8), P)    # x > h -> A
    ok = True
    ok &= check("scalar_step branches", (b1, b2) == ("B", "A"))
    ok &= check("scalar_step B exact", x1 == F(1, 10) + F(1, 4) *
                (P["c"] - F(1, 10)) * P["q"])
    ok &= check("scalar_step A exact", x2 == F(3, 4) * F(1, 8) +
                F(1, 4) * P["L"])
    return ok


def t_word_detect():
    # synthetic: (AAB)^inf; the string is periodic from position 1 as (BAA)^inf
    it = list("AB" + "AAB" * 12)
    i0, n, w = detect_eventual_word(it, rep=10)
    ok = True
    ok &= check("detect_eventual_word i0", i0 == 1)
    ok &= check("detect_eventual_word word", w == "BAA")
    # truncated last run must not break detection
    it2 = list("AB" + "AAB" * 12 + "AA")
    i0b, nb, wb = detect_eventual_word(it2, rep=10)
    ok &= check("detect_eventual_word truncated tail", (i0b, wb) == (1, "BAA"))
    return ok


def t_certificates():
    ok = True
    # K=9 canonical: word B A^11, certified
    P = boundary_params(9)
    p = replay.exact_bin_probs(replay.near_cdf_exact, 9)
    g = replay.exact_bin_probs(replay.xo_cdf_exact, 9)
    for t in range(10):
        D = replay.cumsum_excess(p, 9)
        s = replay.argmax_first(D)
        if s in (7, 8):
            break
        p = replay.update_state(p, g, s, ALPHA)
    x0 = p[7]
    xs, it = scalar_orbit(x0, P, 600)
    i0, n, w = detect_eventual_word(it, rep=20)
    ok &= check("K=9 word is B A^11", w == "B" + "A" * 11)
    c = certify_word(9, x0, w, xs=xs, it=it, i0=i0)
    ok &= check("K=9 certificate ok", bool(c.get("ok")))
    ok &= check("K=9 certificate entry <= 40",
                c.get("t_entry") is not None and c["t_entry"] <= 40)
    # committed census row agreement for K=105
    rows = list(csv.DictReader(
        (REPO / "experiments" / "v03_t7_cyclicity" / "results" /
         "t7_census_certificates.csv").open(encoding="utf-8")))
    row105 = next(r for r in rows if r["K"] == "105")
    P105 = boundary_params(105)
    p = replay.exact_bin_probs(replay.near_cdf_exact, 105)
    g = replay.exact_bin_probs(replay.xo_cdf_exact, 105)
    for t in range(10):
        D = replay.cumsum_excess(p, 105)
        s = replay.argmax_first(D)
        if s in (P105["m"] - 1, P105["m"]):
            break
        p = replay.update_state(p, g, s, ALPHA)
    xs, it = scalar_orbit(p[P105["m"] - 1], P105, 800)
    i0, n, w = detect_eventual_word(it, rep=20)
    ok &= check("K=105 word matches committed census", w == row105["word"])
    ok &= check("K=105 period 19", n == 19)
    return ok


def t_windows():
    ok = True
    ok &= check("pure window K=9 -> 11", pure_window_membership(9) == 11)
    ok &= check("pure window K=15 -> 10", pure_window_membership(15) == 10)
    ok &= check("pure window K=21 -> 9", pure_window_membership(21) == 9)
    ok &= check("pure window K=93 -> 9", pure_window_membership(93) == 9)
    ok &= check("pure window K=99 mixed", pure_window_membership(99) is None)
    ok &= check("pure window K=189 mixed",
                pure_window_membership(189) is None)
    ok &= check("pure window K=195 -> 8", pure_window_membership(195) == 8)
    ok &= check("pure window K=3003 -> 8", pure_window_membership(3003) == 8)
    ok &= check("sigma_3 monotone probe",
                sigma_strictly_monotone(9, 1203, 3))
    ok &= check("sigma_4 monotone probe",
                sigma_strictly_monotone(10, 1204, 4))
    return ok


def t_alt4():
    ok = True
    for K in (16, 100, 598, 1000, 4000):
        P = boundary_params(K)
        beta, s = P["beta"], P["sigma"]
        thr = (F(4, 3) - s) / beta
        m1 = F(3, 4) * (beta * F(3, 4) + s) - thr
        m2 = thr - F(9, 16) * (beta + s)
        ok &= check(f"r=4 alternation K={K}", m1 > 0 and m2 >= 0)
    return ok


def t_lock():
    ok = True
    for K in (7, 12, 13, 14):
        P = boundary_params(K)
        m, a, L = P["m"], P["a"], P["L"]
        p = replay.exact_bin_probs(replay.near_cdf_exact, K)
        g = replay.exact_bin_probs(replay.xo_cdf_exact, K)
        lastB = None
        t_entry = None
        good = True
        gap_prev = None
        for t in range(60):
            D = replay.cumsum_excess(p, K)
            s = replay.argmax_first(D)
            if s in (m - 1, m) and t_entry is None:
                t_entry = t
            if s == m - 1:
                lastB = t
            t_lock = lastB if lastB is not None else t_entry
            if t_lock is not None and t > t_lock:
                if s != m or D[s - 1] != a:
                    good = False
                gap = abs(p[m - 1] - L)
                if gap_prev is not None and gap > gap_prev:
                    good = False
                gap_prev = gap
            p = replay.update_state(p, g, s, ALPHA)
        ok &= check(f"LOCK K={K} sel=m & d=a & |x-L| non-increasing "
                    f"after lock (covers never-B case)", good)
    return ok


def t_theorem_b():
    import run_b
    ok = True
    ok &= check("Theta(9)=41/45", run_b.theta(9) == F(41, 45))
    ok &= check("Theta(6)=43/42", run_b.theta(6) == F(43, 42))
    p0 = run_b.make_init(21, F(9, 10))
    P = boundary_params(21)
    m = P["m"]
    C = sum(p0[:m])
    # A1 verified inside make_init; check conservation setup
    ok &= check("B make_init F_{0,m}=C", C == F(9, 10))
    r = run_b.run_case(21, F(9, 10), 400)
    ok &= check("B case invariants (K=21,C=0.9)",
                r["sel_le_m"] and r["conserved"] and r["scalar_ok"]
                and r["both_branches"])
    r2 = run_b.run_case(12, F(9, 10), 400)
    ok &= check("B case LOCK (K=12,C=0.9)",
                r2["sel_le_m"] and r2["eventual_m"])
    return ok


def t_theorem_c():
    import run_c
    ok = True
    b = run_c.bounds(9, F(1, 4))
    ok &= check("C bounds K=9 alpha=1/4 NA=11", b["NA"] == 11)
    b2 = run_c.bounds(9, F(1, 8))
    ok &= check("C bounds K=9 alpha=1/8 NA=18", b2["NA"] == 18)
    ok &= check("C bounds K=9 alpha=1/8 NA>=K (freeze-rule risk)",
                b2["NA"] >= 9)
    # counterexample reproduction (short)
    t_star, s_list = run_c.entry_and_itinerary(9, F(1, 8), 90)
    m = run_c.params(9)[0]
    fr = run_c.classify_frozen(s_list, 9, 90)
    ok &= check("C counterexample frozen=LOCK", fr == "LOCK")
    tau0 = t_star + max(b2["NB"], b2["U1"])
    two = "CYCLE" if any(s == m - 1 for s in s_list[tau0 + 1:]) else "LOCK"
    ok &= check("C counterexample two-test=CYCLE", two == "CYCLE")
    return ok


def t_family_cover():
    """Regression: the committed fixed family-cover certificate verifies
    by fixed-parameter re-computation with the envelope constant bound to
    the polynomial certificate; tampered certificates are REJECTED
    (negative tests, review R3 directive)."""
    import copy
    import importlib.util
    import json
    cert_path = (REPO / "experiments" / "v03_t7_cyclicity" / "results" /
                 "t7_family_cover_cert.json")
    if not cert_path.exists():
        return check("family cover cert present", False)
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    spec = importlib.util.spec_from_file_location(
        "run_t7", str(REPO / "experiments" / "v03_t7_cyclicity" / "run_t7.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    res = m.verify_family_cover()
    ok = True
    ok &= check("family cover certificate verified (fixed-parameter "
                "re-computation)", res["verified"])
    ok &= check("family cover K0 committed value", res["K0"] == cert["K0"])
    ok &= check("family cover strip positive",
                F(cert["th_hi"]) > F(cert["th_lo"]))
    ok &= check("family cover island has 5 phases", len(cert["island"]) == 5)
    ok &= check("family cover sequence W5", cert["sequence"] == [2, 2, 3, 2, 2])
    ok &= check("family cover margins positive",
                F(cert["min_inclusion_margin"]) > 0 and
                F(cert["min_crossing_margin"]) > 0)

    def rejected(name, mutated):
        r = m.verify_family_cover(cert=copy.deepcopy(mutated), hard=False)
        return check(f"family cover REJECTS {name}",
                     r is None or r.get("verified") is False)

    # negative tests: wrong envelope (self-consistent-but-fake C5), wrong
    # box fields, wrong island, wrong sequence must all be rejected
    fake = copy.deepcopy(cert)
    fake["C5_per_n"] = "3056989494373775/1712558184037213392"
    ok &= rejected("wrong envelope constant (fake small C5)", fake)
    for field, value in (("beta_lo", "-1"), ("sig_hi", "99"),
                         ("band_hi", "-1"), ("dmin_limit", "0"),
                         ("n_pop", 0), ("strip_width", "0")):
        bad = copy.deepcopy(cert)
        bad[field] = value
        ok &= rejected(f"tampered {field}", bad)
    bad_seq = copy.deepcopy(cert)
    bad_seq["sequence"] = [3, 3, 3, 3, 3]
    ok &= rejected("wrong sequence", bad_seq)
    bad_isl = copy.deepcopy(cert)
    bad_isl["island"] = [["0", "1/1000"]] * 5
    ok &= rejected("wrong island", bad_isl)
    bad_sigma_inf = copy.deepcopy(cert)
    bad_sigma_inf["sigma_inf"] = "0"
    ok &= rejected("wrong sigma_inf", bad_sigma_inf)
    return ok


def t_theorem_d():
    """Regression for the stopping-rule experiment (reviews R4/R5):
    real, failure-capable checks only — Lambda recomputed from the word
    string, periodic-point branch back-substitution, per-period increment
    ratio consistency on the actual trajectory, K=6 rule exception, LOCK
    identity. No tautologies."""
    import run_d
    ok = True
    for K in (9, 10, 11):
        r = run_d.analyze_cycle(K)
        ok &= check(f"D CYCLE K={K}: best_t=2 + tail cert",
                    r["best_t"] == 2 and r["cycle_cert"])
        # Lambda recomputed from the WORD STRING: (3/4)^#A * beta^#B
        P = boundary_params(K)
        lam_direct = (F(3, 4) ** r["word"].count("A")) * \
                     (P["beta"] ** r["word"].count("B"))
        ok &= check(f"D CYCLE K={K}: Lambda = (3/4)^#A(w) beta^#B(w) "
                    f"recomputed from word", F(r["lam_word"]) == lam_direct)
        # periodic point realizes the word branch-by-branch
        xstar = F(r["xstar"])
        x = xstar
        side_ok = True
        for br in r["word"]:
            side_ok &= ((x > P["h"]) == (br == "A"))
            x = run_d_forced(x, br, P)
        ok &= check(f"D CYCLE K={K}: fixed point realizes word + returns",
                    side_ok and x == xstar)
        # per-period increment ratio on the ACTUAL trajectory at the
        # trajectory's own phase at t0 (found by matching selections):
        # x_{t0+n} - ps[phi] == Lambda * (x_{t0} - ps[phi]), exactly
        states, t_star, P2, g = run_d.trajectory(K, 300)
        m = P["m"]
        t0 = len(states) - 1 - r["n"]
        sel_at = lambda t: replay.argmax_first(
            replay.cumsum_excess(states[t], K))
        span = range(len(states) - t0)
        phi = next(cand for cand in range(r["n"])
                   if all(sel_at(t0 + i) ==
                          (m if r["word"][(cand + i) % r["n"]] == "A"
                           else m - 1)
                          for i in span))
        x_phi = xstar
        for br in r["word"][:phi]:
            x_phi = run_d_forced(x_phi, br, P)
        lhs = states[t0 + r["n"]][m - 1] - x_phi
        rhs = F(r["lam_word"]) * (states[t0][m - 1] - x_phi)
        ok &= check(f"D CYCLE K={K}: x-increment ratio == Lambda exactly "
                    f"(own phase)", lhs == rhs)
    # K=9 word Lambda exact value (reviewer's independent computation)
    r9 = run_d.analyze_cycle(9)
    ok &= check("D K=9 Lambda = 688905/16777216 (reviewer value)",
                F(r9["lam_word"]) == F(688905, 16777216))
    # LOCK identity spot: K=7
    r7 = run_d.analyze_lock(7)
    ok &= check("D LOCK K=7: best_t=2 + quad cert + S'<0",
                r7["best_t"] == 2 and r7["quad_cert"] and r7["S_sign"] == "neg")
    # K=6 online-rule exception: stop-at-first-boundary is NOT optimal
    states, t_star, P, g = run_d.trajectory(6, 60)
    u1 = run_d.U(states[1], 6)
    u2 = run_d.U(states[2], 6)
    ok &= check("D K=6: t*=1 but U(p1) > U(p2) (rule exception)",
                t_star == 1 and u1 > u2 and u2 == F(335, 2304))
    return ok


def run_d_forced(x, br, P):
    if br == "A":
        return F(3, 4) * x + P["L"] * F(1, 4)
    return x + (P["c"] - x) * P["q"] * F(1, 4)


def t_a_plateau_allk():
    """Research A regression: the committed all-K plateau polynomial
    certificates re-verify from JSON (fixed-parameter re-computation +
    wp1_rows cross-check + zone-set completeness/uniqueness + source-hash
    binding); tampered certificates must be REJECTED."""
    import importlib.util
    import json
    import tempfile
    spec = importlib.util.spec_from_file_location(
        "run_a", str(REPO / "experiments" / "v03_p_plateau_allk" / "run_a.py"))
    run_a = importlib.util.module_from_spec(spec)
    sys.modules["run_a"] = run_a        # for later attribute imports
    spec.loader.exec_module(run_a)
    ok = True
    fails = run_a.verify_from_json()
    ok &= check("A all-K plateau certificates verified from JSON "
                "(3 zones, 12 polynomials, zone set exact)", fails == [])
    cert_path = (REPO / "experiments" / "v03_p_plateau_allk" / "results" /
                 "a_plateau_cert.json")
    if not cert_path.exists():
        return check("A certificate file present", False)
    data = json.loads(cert_path.read_text(encoding="utf-8"))

    def verify_with(data_obj):
        with tempfile.TemporaryDirectory() as td:
            (Path(td) / "a_plateau_cert.json").write_text(
                json.dumps(data_obj), encoding="utf-8")
            old = run_a.RESULTS
            run_a.RESULTS = Path(td)
            try:
                return run_a.verify_from_json()
            finally:
                run_a.RESULTS = old

    bad = json.loads(json.dumps(data))
    co = bad["zones"][0]["certificates"]["right_num"]["shift_coeffs"]
    co[0] = "-" + co[0] if not co[0].startswith("-") else co[0][1:]
    ok &= check("A verifier REJECTS tampered shift coefficient",
                verify_with(bad) != [])
    ok &= check("A verifier REJECTS empty zone set",
                verify_with({"zones": []}) != [])
    bad2 = json.loads(json.dumps(data))
    bad2["zones"] = [z for z in bad2["zones"] if z["r"] != 5]
    ok &= check("A verifier REJECTS missing r=5 zone",
                verify_with(bad2) != [])
    bad3 = json.loads(json.dumps(data))
    bad3["source_sha256"] = "0" * 64
    ok &= check("A verifier REJECTS source-hash mismatch",
                verify_with(bad3) != [])
    bad4 = json.loads(json.dumps(data))
    z0 = bad4["zones"][0]
    z0["certificates"]["left_num"]["finite"] = \
        z0["certificates"]["left_num"]["finite"][:-1]
    if z0["certificates"]["left_num"]["finite"]:
        ok &= check("A verifier REJECTS truncated finite prefix",
                    verify_with(bad4) != [])
    else:
        # zone starts at N=n0 (empty prefix is legitimate); fabricate a
        # gap instead by dropping the start binding
        bad5 = json.loads(json.dumps(data))
        for z in bad5["zones"]:
            for cname in z["certificates"]:
                z["certificates"][cname]["covers_from"] = \
                    z["certificates"][cname]["N"] + 1
        ok &= check("A verifier REJECTS start-coverage gap",
                    verify_with(bad5) != [])
    # owner review 8: SELF-CONSISTENT start shift -- move covers_from AND
    # N together by +100, recompute the shift coefficients at the new N,
    # empty the prefix.  Internally consistent, but the first 100 zone
    # parameters are silently uncovered: must be rejected because the
    # REQUIRED start (zone n0) comes from the caller, not the certificate.
    _shift = run_a.pp_shift
    bad6 = json.loads(json.dumps(data))
    for z in bad6["zones"]:
        r = z["r"]
        p = int(z["rho"].split("/")[0])
        qd = int(z["rho"].split("/")[1])
        _, left, right = run_a.comparisons_rf(r, p, qd)
        polys = {"left_num": left.num, "left_den": left.den,
                 "right_num": right.num, "right_den": right.den}
        for cname, P in polys.items():
            c = z["certificates"][cname]
            c["N"] = c["N"] + 100
            c["covers_from"] = c["N"]
            c["finite"] = []
            c["shift_coeffs"] = [str(x) for x in _shift(P, c["N"])]
    ok &= check("A verifier REJECTS self-consistent start shift "
                "(first 100 zone values uncovered)", verify_with(bad6) != [])
    # zone sanity: the three infinite zones match the committed table
    zones = [(z["r"], z["rho"], z["n0"]) for z in data["zones"]]
    ok &= check("A zones are (3,1/9,32),(4,2/5,2),(5,11/16,10)",
                zones == [(3, "1/9", 32), (4, "2/5", 2), (5, "11/16", 10)])
    return ok


def t_b_lock_allk():
    """Research B regression (review-7 fixed): committed certificates for
    ALL THREE LOCK residues re-verify from JSON with set completeness;
    tampering must be rejected; the K=12 owner anchors and the
    second-selection closed form hold."""
    import importlib.util
    import json
    import tempfile
    spec = importlib.util.spec_from_file_location(
        "run_b2", str(REPO / "experiments" / "v03_b_lock_allk" / "run_b2.py"))
    run_b2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run_b2)
    ok = True
    fails = run_b2.verify_from_json()
    ok &= check("B LOCK certificates verified from JSON "
                "(96 sub-families x 6-7 claims, set exact)", fails == [])
    cert_path = (REPO / "experiments" / "v03_b_lock_allk" / "results" /
                 "b_lock_cert.json")
    if not cert_path.exists():
        return check("B certificate file present", False)
    data = json.loads(cert_path.read_text(encoding="utf-8"))
    ok &= check("B 96 sub-families stored (3 residues x 32)",
                len(data["families"]) == 96 and
                len({(f["r"], f["nu"]) for f in data["families"]}) == 96)
    ok &= check("B old wrong-direction claim C5 absent, repair C5p present",
                all("C5p_-(Sp+Ep)>0" in f["claims"] and
                    "C5_tail>=q(star)" not in f["claims"]
                    for f in data["families"]))
    ok &= check("B crossing claims present (P1 L>h; P2 x2/P3 x3 by residue)",
                all("P1_L-h>0" in f["claims"] and
                    ("P2_x2-h>0" in f["claims"]) == (f["r"] == 0) and
                    ("P2_h-x2>0" in f["claims"] and "P3_x3-h>0" in
                     f["claims"]) == (f["r"] in (1, 2))
                    for f in data["families"]))
    # owner review 8: all-K boundary-entry lemma (analytic K >= 33)
    ok &= check("B lemma t*=2 verified (K>=33 analytic)", run_b2.lemma_t2()
                == [])
    # owner review 8: a certificate valid only from u=3 must be REJECTED
    # when the theorem requires the claim from u=0 (P(u) = u - 2 < 0 there)
    _ra = sys.modules["run_a"]          # imported by run_b2 during its load
    _RF = _ra.RF
    _shift = _ra.pp_shift
    rf_u2 = _RF([F(-2), F(1)], [F(1)])          # u - 2
    cert_u2 = {"strict": True, "N": 3, "covers_from": 3,
               "shift_coeffs": [str(x) for x in _shift(rf_u2.num, 3)],
               "finite": []}
    stored_u2 = {"num": cert_u2, "den_sign": 1,
                 "den_cert": {"strict": True, "N": 3, "covers_from": 3,
                              "shift_coeffs": ["1"], "finite": []}}
    ok &= check("B verify_claim REJECTS u-2 cert when start u=0 required",
                run_b2.verify_claim(rf_u2, stored_u2, 0, True) != [])
    ok &= check("B verify_claim ACCEPTS same cert for matching start u=3",
                run_b2.verify_claim(rf_u2, stored_u2, 3, True) == [])

    def verify_with(data_obj):
        with tempfile.TemporaryDirectory() as td:
            (Path(td) / "b_lock_cert.json").write_text(
                json.dumps(data_obj), encoding="utf-8")
            old = run_b2.RESULTS
            run_b2.RESULTS = Path(td)
            try:
                return run_b2.verify_from_json()
            finally:
                run_b2.RESULTS = old

    bad = json.loads(json.dumps(data))
    co = bad["families"][5]["claims"]["C1_-Sp>0"]["num"]["shift_coeffs"]
    co[0] = "-" + co[0] if not co[0].startswith("-") else co[0][1:]
    ok &= check("B verifier REJECTS tampered shift coefficient",
                verify_with(bad) != [])
    bad2 = json.loads(json.dumps(data))
    nu0 = next(f for f in bad2["families"] if f["r"] == 0 and f["nu"] == 0)
    bad2["families"] = [dict(nu0, nu=f["nu"], u0=f["u0"])
                        if f["r"] == 0 else f for f in bad2["families"]]
    ok &= check("B verifier REJECTS 32 clones of nu=0 replacing r=0 set",
                verify_with(bad2) != [])
    bad3 = json.loads(json.dumps(data))
    bad3["families"] = bad3["families"][:50]
    ok &= check("B verifier REJECTS missing sub-families",
                verify_with(bad3) != [])
    # owner review 8: self-consistent start shift on ONE family's claims
    # (covers_from and N moved together, coefficients recomputed by the
    # verifier's own polynomials, prefix emptied) must be rejected
    _shift = sys.modules["run_a"].pp_shift
    bad4 = json.loads(json.dumps(data))
    entry = next(f for f in bad4["families"] if (f["r"], f["nu"]) == (0, 5))
    sub = next(s for s in run_b2.subfamilies()
               if (s["r"], s["nu"]) == (0, 5))
    Dsub = run_b2.build_family(sub)
    for name, (rf_, strict) in run_b2.claims_for(Dsub).items():
        c = entry["claims"][name]["num"]
        c["N"] = c["N"] + 5
        c["covers_from"] = c["N"]
        c["finite"] = []
        num = rf_.num if entry["claims"][name]["den_sign"] == 1 else \
            [-x for x in rf_.num]
        c["shift_coeffs"] = [str(x) for x in _shift(num, c["N"])]
        cd = entry["claims"][name]["den_cert"]
        cd["N"] = cd["N"] + 5
        cd["covers_from"] = cd["N"]
        cd["finite"] = []
        den = rf_.den if entry["claims"][name]["den_sign"] == 1 else \
            rf_.den * (-1)
        cd["shift_coeffs"] = [str(x) for x in _shift(den, cd["N"])]
    ok &= check("B verifier REJECTS self-consistent start shift "
                "(first 5 sub-family values uncovered)",
                verify_with(bad4) != [])
    # owner anchors at K=12: old C5 false, vertex > 1, repair positive
    D12 = run_b2.build_family(next(s for s in run_b2.subfamilies()
                                   if s["r"] == 0 and s["nu"] == 2))
    _pe = sys.modules["run_a"].pp_eval
    u = 0
    Sp = _pe(D12["Sp"].num, u) / _pe(D12["Sp"].den, u)
    Ep = _pe(D12["Ep"].num, u) / _pe(D12["Ep"].den, u)
    U2 = _pe(D12["U2"].num, u) / _pe(D12["U2"].den, u)
    Ui = _pe(D12["Uinf"].num, u) / _pe(D12["Uinf"].den, u)
    ok &= check("B K=12 owner anchor: old C5 value > 0 (false claim)",
                Sp * Sp + (U2 - Ui) * Ep == F(205209, 1073741824))
    ok &= check("B K=12 owner anchor: vertex -S'/E' > 1 (outside range)",
                -Sp / Ep == F(226240, 157449))
    ok &= check("B K=12 owner anchor: -(S'+E') > 0 (repair claim)",
                -(Sp + Ep) > 0)
    # second-selection closed form on a replay sample (all three residues,
    # tie points included: even-K tie K=96; odd K examples)
    import replay as rp
    for K in (12, 24, 96, 25, 31, 61, 8, 14, 68, 194):
        p0 = rp.exact_bin_probs(rp.near_cdf_exact, K)
        g = rp.exact_bin_probs(rp.xo_cdf_exact, K)
        n, r = divmod(K, 6)
        s0 = 3 * n + (1 if r == 2 else 0)
        p = rp.update_state(p0, g, s0, F(1, 4))
        s1 = rp.argmax_first(rp.cumsum_excess(p, K))
        ok &= check(f"B s1 closed form K={K}", s1 == run_b2.s1_closed(K))
    return ok


def t_c_alpha_endpoints():
    """Research C regression (review 9): the SINGLE complete acceptance
    function run_c2.verify_table is used for the positive check and for
    every single-break negative (missing record, duplicate record,
    deleted brackets, wrong bracket inside legal range, wrong degree,
    wrong hash); Sturm counterexample self-test and the G-ep anchor
    are re-derived exactly."""
    import importlib.util
    import json
    spec = importlib.util.spec_from_file_location(
        "run_c2", str(REPO / "experiments" / "v03_c_alpha_endpoints" /
                      "run_c2.py"))
    run_c2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run_c2)
    ok = True
    ok &= check("C Sturm self-test (x^2(x-1/2) -> 1 on (0,1) etc.)",
                run_c2.sturm_selftest())
    res = run_c2.classify(10, 2, 5)
    ok &= check("C anchor K=10 rho=2/5: N^+ == 1",
                res["plus"]["roots_open_0_1"] == 1)
    ok &= check("C anchor K=10 rho=2/5: N^- == 1",
                res["minus"]["roots_open_0_1"] == 1)
    lo = F(res["plus"]["brackets"][0]["bracket"][0])
    hi = F(res["plus"]["brackets"][0]["bracket"][1])
    ok &= check("C anchor bracket inside G-ep 0.245996/0.245997",
                F(245996, 10 ** 6) < lo and hi < F(245997, 10 ** 6))
    data = json.loads((REPO / "experiments" / "v03_c_alpha_endpoints" /
                       "results" / "c_alpha_endpoints.json").read_text(
                           encoding="utf-8"))
    ok &= check("C committed table PASSES the single complete acceptance "
                "(verify_table: set/counts/degrees/brackets/isolation/"
                "genuine/sha)", run_c2.verify_table(data) == [])

    def tampered(**mut):
        bad = json.loads(json.dumps(data))
        mut(bad)
        return bad

    def drop_record(d):
        d["targets"].remove(d["targets"][3])

    def dup_record(d):
        d["targets"].append(dict(d["targets"][0]))

    def del_brackets(d):
        for t in d["targets"]:
            t["plus"]["brackets"] = []

    def wrong_bracket(d):
        # bracket inside the legal range (0,1) but isolating NO root --
        # every other field untouched
        t = next(x for x in d["targets"] if x["plus"]["brackets"])
        t["plus"]["brackets"][0]["bracket"] = ["1/100", "1/50"]

    def wrong_degree(d):
        d["targets"][0]["plus"]["num_degree"] += 1

    def wrong_hash(d):
        d["source_sha256"] = "0" * 64

    for name, mut in (("missing record", drop_record),
                      ("duplicate record", dup_record),
                      ("deleted brackets", del_brackets),
                      ("wrong bracket in legal range", wrong_bracket),
                      ("wrong degree", wrong_degree),
                      ("wrong source hash", wrong_hash)):
        bad = json.loads(json.dumps(data))
        mut(bad)
        ok &= check(f"C acceptance REJECTS {name} (single break)",
                    run_c2.verify_table(bad) != [])
    return ok


def t_e_cycle_allk():
    """Theorem E regression: committed CYCLE certificates re-verify from
    JSON with exact set/claim binding; single-break negatives (missing
    residue family, dropped K=10 exception record, self-consistent start
    shift, strictness downgrade, tampered key polynomial, K=12-style
    anchors); independent second implementation on a small sample."""
    import importlib.util
    import json
    import tempfile
    spec = importlib.util.spec_from_file_location(
        "run_e", str(REPO / "experiments" / "v03_e_cycle_allk" / "run_e.py"))
    run_e = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run_e)
    sys.modules.setdefault("run_e", run_e)
    ok = True
    ok &= check("E lemma t*=2 (CYCLE, anchor D_{2,m-1}, M(33) bound) "
                "verified", run_e.lemma_t2_cycle() == [])
    fails = run_e.verify_from_json()
    ok &= check("E CYCLE certificates verified from JSON "
                "(96 sub-families x 6 square-completion claims + "
                "finite set via D, set exact)", fails == [])
    cert_path = (REPO / "experiments" / "v03_e_cycle_allk" / "results" /
                 "e_cycle_cert.json")
    if not cert_path.exists():
        return check("E certificate file present", False)
    data = json.loads(cert_path.read_text(encoding="utf-8"))
    ok &= check("E 96 sub-families stored (3 residues x 32)",
                len(data["families"]) == 96 and
                len({(f["r"], f["nu"]) for f in data["families"]}) == 96)
    ok &= check("E config records the review-10 identity and route",
                data.get("config", {}).get("identity") ==
                "U = Q + 2 z Phi + K E2 z^2" and
                data.get("config", {}).get("route") ==
                "square completion (S1-S4)")
    ok &= check("E finite exception set is exactly "
                "{9,10,11,16,17,23} via D certificates",
                [e["K"] for e in data.get("finite_exceptions", [])] ==
                [9, 10, 11, 16, 17, 23] and
                all("verify_d_independent" in e["via"]
                    for e in data["finite_exceptions"]))
    ok &= check("E claim set is the 6 square-completion claims "
                "(no legacy 5-point Phi/R claims)",
                all(set(f["claims"]) ==
                    {"C0a_U0-U2>=0", "C0b_U1-U2>=0", "S1_A>0",
                     "S2_H2>0", "S3_-Hp(3/4)>0", "S4_H(3/4)>=0"}
                    for f in data["families"]))

    def verify_with(data_obj):
        with tempfile.TemporaryDirectory() as td:
            (Path(td) / "e_cycle_cert.json").write_text(
                json.dumps(data_obj), encoding="utf-8")
            old = run_e.RESULTS
            run_e.RESULTS = Path(td)
            try:
                return run_e.verify_from_json()
            finally:
                run_e.RESULTS = old

    # (a) missing residue family r=5
    bad = json.loads(json.dumps(data))
    bad["families"] = [f for f in bad["families"] if f["r"] != 5]
    ok &= check("E verifier REJECTS missing residue family r=5",
                verify_with(bad) != [])
    # (b) dropped finite-exception record
    bad = json.loads(json.dumps(data))
    bad["finite_exceptions"] = bad["finite_exceptions"][:-1]
    ok &= check("E verifier REJECTS dropped finite-set record",
                verify_with(bad) != [])
    # (b2) finite record without the D reference
    bad = json.loads(json.dumps(data))
    bad["finite_exceptions"][0]["via"] = "600-step window check"
    ok &= check("E verifier REJECTS finite record lacking D reference",
                verify_with(bad) != [])
    # (c) self-consistent start shift on one family
    _shift = sys.modules["run_a"].pp_shift
    bad = json.loads(json.dumps(data))
    entry = next(f for f in bad["families"] if (f["r"], f["nu"]) == (3, 5))
    sub = next(s for s in run_e.subfamilies()
               if (s["r"], s["nu"]) == (3, 5))
    Dsub = run_e.build_family(sub)
    for name, (rf_, strict) in run_e.claims_for(Dsub).items():
        c = entry["claims"][name]["num"]
        c["N"] += 5
        c["covers_from"] = c["N"]
        c["finite"] = []
        num = rf_.num if entry["claims"][name]["den_sign"] == 1 else \
            [-x for x in rf_.num]
        c["shift_coeffs"] = [str(x) for x in _shift(num, c["N"])]
        cd = entry["claims"][name]["den_cert"]
        cd["N"] += 5
        cd["covers_from"] = cd["N"]
        cd["finite"] = []
        den = rf_.den if entry["claims"][name]["den_sign"] == 1 else \
            rf_.den * (-1)
        cd["shift_coeffs"] = [str(x) for x in _shift(den, cd["N"])]
    ok &= check("E verifier REJECTS self-consistent start shift",
                verify_with(bad) != [])
    # (d) strictness downgrade on a strict claim
    bad = json.loads(json.dumps(data))
    entry = next(f for f in bad["families"] if (f["r"], f["nu"]) == (4, 3))
    c = entry["claims"]["S1_A>0"]["num"]
    c["strict"] = False
    ok &= check("E verifier REJECTS strictness downgrade",
                verify_with(bad) != [])
    # (e) tampered key coefficient (H(3/4) claim)
    bad = json.loads(json.dumps(data))
    entry = next(f for f in bad["families"] if (f["r"], f["nu"]) == (5, 7))
    co = entry["claims"]["S4_H(3/4)>=0"]["num"]["shift_coeffs"]
    co[0] = "-" + co[0] if not co[0].startswith("-") else co[0][1:]
    ok &= check("E verifier REJECTS tampered S4 coefficient",
                verify_with(bad) != [])
    # (f) missing transient claim (drop C0a everywhere)
    bad = json.loads(json.dumps(data))
    for f in bad["families"]:
        del f["claims"]["C0a_U0-U2>=0"]
    ok &= check("E verifier REJECTS missing transient claim C0a",
                verify_with(bad) != [])
    # (g) legacy route marker (old false identity) must be rejected
    bad = json.loads(json.dumps(data))
    bad["config"]["identity"] = "U = Q + z Phi + K E2 z^2"
    ok &= check("E verifier REJECTS the legacy factor-1 identity marker",
                verify_with(bad) != [])
    # independent second implementation on a small exact sample
    spec2 = importlib.util.spec_from_file_location(
        "vei", str(REPO / "experiments" / "v03_e_cycle_allk" /
                   "verify_e_independent.py"))
    vei = importlib.util.module_from_spec(spec2)
    spec2.loader.exec_module(vei)
    bad_k = [res for res in (vei.analyze(K) for K in
                             (9, 10, 11, 15, 16, 17, 23, 63, 124, 191,
                              195, 196, 197, 999, 1000, 1001)) if res]
    ok &= check("E independent engine (review-10 repaired): theorem + "
                "decomp + IDENTITY BACK-SUBSTITUTION + S1-S4 on 16-K "
                "sample (own model definitions)", not bad_k)
    if bad_k:
        print("   independent failures:", bad_k)
    # Recompute both owner counterexamples from the independent model.
    # Phi is the half cross-term, so the old factor-1 identity and old
    # lower bound must reproduce the exact negative discrepancies.
    K9 = 9
    p9 = [vei.bins_p0(K9)]
    g9 = vei.bins_g(K9)
    for _ in range(3):
        p9.append(vei.step(p9[-1], g9, K9)[0])
    P9 = boundary_params(K9)
    m9, h9 = P9["m"], F(1, K9)
    gamma9 = lambda x: (P9["c"] - x) / (1 - P9["q"])
    x29, x39 = p9[2][m9 - 1], p9[3][m9 - 1]
    e9 = [p9[2][i] - gamma9(x29) * g9[i] for i in range(m9 - 1)]

    def Q9(x):
        profile = ([gamma9(x) * gi for gi in g9[:m9 - 1]]
                   + [x] + p9[0][m9:])
        return K9 * sum((pi - h9) ** 2 for pi in profile)

    def Phi9(x):
        return K9 * sum((gamma9(x) * g9[i] - h9) * e9[i]
                        for i in range(m9 - 1))

    E9 = K9 * sum(ei * ei for ei in e9)
    U29, U39 = [K9 * sum((pi - h9) ** 2 for pi in p9[t])
                for t in (2, 3)]
    old_identity_gap = U29 - (Q9(x29) + Phi9(x29) + E9)
    old_R39 = Q9(x39) + F(3, 4) * Phi9(x39) - U29
    old_bound_gap = (U39 - U29) - old_R39
    ok &= check("E owner anchors: independent K=9 replay reproduces "
                "both old-formula counterexamples exactly",
                old_identity_gap == -F(2459287, 10077696) and
                old_bound_gap == -F(12642529, 241864704))
    return ok


def t_bii_step_bound():
    """Theorem B-II last-two-cells step bound (closeout item 1, 2026-10-19):
    h/eps < 4/(2(b-h)-1) = 12K/(2K-1) <= 164/27 < 7 for r=5, K>=41,
    where the middle identity uses 2m-2-K = (2K-1)/3. The old intermediate
    inequality 4/(2(6/7-1/41)-1) had the wrong direction (2(C-h)-1 is
    increasing in C; the bound must use C > b, not C <= 6/7)."""
    ok = True
    for K in (41, 47, 53, 65, 101, 479, 1001):
        P = boundary_params(K)
        m, b, h = P["m"], P["b"], P["h"]
        ok &= check(f"B-II r=5 identity 2m-2-K=(2K-1)/3 K={K}",
                    P["r"] == 5 and 2 * m - 2 - K == (2 * K - 1) // 3)
        # 4/(2(b-h)-1) == 12K/(2K-1) exactly
        den = 2 * (b - h) - 1
        ok &= check(f"B-II 4/(2(b-h)-1) = 12K/(2K-1) K={K}",
                    F(4) / den == F(12 * K, 2 * K - 1))
    # monotone decreasing bound, max at K=41, uniform 7-step conclusion
    ok &= check("B-II bound form 12K/(2K-1) = 6 + 6/(2K-1)",
                all(F(12 * k, 2 * k - 1) == 6 + F(6, 2 * k - 1)
                    for k in (41, 47, 53, 1001)))
    ok &= check("B-II max at K=41 equals 164/27",
                F(12 * 41, 2 * 41 - 1) == F(164, 27))
    ok &= check("B-II 164/27 < 7 hence floor(h/eps)+1 <= 7",
                F(164, 27) < 7)
    # per-K upper bounds all give floor+1 = 7 (corrects old '7/6/6' record
    # which was computed from the wrong-direction inequality)
    ok &= check("B-II K=41/53/65 bounds 164/27, 212/35, 260/43 all > 6",
                F(164, 27) > 6 and F(212, 35) > 6 and F(260, 43) > 6)
    return ok


def t_wp4_error_semantics():
    """WP4 three-way status semantics (closeout item 2, 2026-10-19): solver
    anomalies are ERROR, never infeasibility or conflict evidence; the
    targeted fault-injection check must pass."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "wp4", str(REPO / "experiments" / "v03_wp4_bopen" / "wp4.py"))
    wp4 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(wp4)
    ok = True
    ok &= check("WP4 fault-injection self-test (fixed-point anomalies "
                "routed to ERROR)", wp4.fault_injection_check())
    # ground truth: lock word feasible at K=7 (r=1), infeasible at K=9 (r=3)
    m7 = wp4.anchor(7)[0]
    ok &= check("WP4 lock word feasible at K=7",
                wp4.solve_word(7, [m7])["status"] == "feasible")
    m9 = wp4.anchor(9)[0]
    ok &= check("WP4 lock word infeasible at K=9",
                wp4.solve_word(9, [m9])["status"] == "infeasible")
    # robust/subsystem classification on an infeasible word: three-way tags
    tag, conf, errs = wp4.robust_and_subsystems(9, (1, m9))
    ok &= check("WP4 robust tag infeasible at K=9 word (1,m)",
                tag == "robust-infeasible")
    ok &= check("WP4 subsystem errors never occur on clean run", errs == [])
    ok &= check("WP4 subsystem conflicts recorded (SEL-only or stronger)",
                len(conf) >= 1)
    # injected anomaly surfaces as robust-error / separate error, not conflicts
    def flaky(solver, fail_on):
        calls = [0]
        def wrapped(K, word, **kw):
            calls[0] += 1
            if calls[0] in fail_on:
                return {"status": "error", "why": "injected"}
            return solver(K, word, **kw)
        return wrapped
    tag_e, conf_e, errs_e = wp4.robust_and_subsystems(
        9, (1, m9), solver=flaky(wp4.solve_word, {1, 5}))
    ok &= check("WP4 injected robust anomaly -> robust-error (not infeasible)",
                tag_e == "robust-error")
    # the errored subsystem mask must appear ONLY in the error list; the
    # conflict list may lose at most that one mask (it turned into an error)
    conf_masks = {c[:3] for c in conf}
    conf_e_masks = {c[:3] for c in conf_e}
    ok &= check("WP4 injected subsystem anomaly -> separate error list, "
                "conflicts lose at most the errored mask",
                len(errs_e) == 1 and errs_e[0][3] == "injected"
                and errs_e[0][:3] not in conf_e_masks
                and conf_masks - conf_e_masks <= {errs_e[0][:3]})
    return ok


def t_d_independent_negative():
    """Negative tests for the completed independent D verifier (R5):
    tampered records must be REJECTED — wrong best_t, wrong tail_lb,
    wrong license_margin, wrong t_star, empty set, missing K."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "vdi", str(REPO / "experiments" / "v03_d_stopping" /
                  "verify_d_independent.py"))
    vdi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vdi)
    import json
    cert_path = (REPO / "experiments" / "v03_d_stopping" / "results" /
                 "d_certificates.json")
    certs = json.loads(cert_path.read_text(encoding="utf-8"))
    ok = True

    def run_modified(label, mutate):
        bad = json.loads(json.dumps(certs))  # deep copy
        mutate(bad)
        # write to the real path? No: patch module's CERTS via tmp file
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".json",
                                         delete=False,
                                         encoding="utf-8") as fh:
            json.dump(bad, fh)
            tmppath = fh.name
        old = vdi.CERTS
        vdi.CERTS = type(old)(tmppath)
        try:
            # capture exit via main() return
            import io
            import contextlib
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = vdi.main()
            rejected = (rc != 0)
        finally:
            vdi.CERTS = old
        return check(f"D verifier REJECTS {label}", rejected)

    def find_cyc(records):
        return next(r for r in records if r["class"] == "CYCLE")

    ok &= run_modified("wrong best_t (0)",
                       lambda b: find_cyc(b).__setitem__("best_t", 0))
    ok &= run_modified("wrong tail_lb (-100)",
                       lambda b: find_cyc(b).__setitem__("tail_lb", "-100"))
    ok &= run_modified(
        "wrong license_margin (999)",
        lambda b: find_cyc(b).__setitem__("license_margin", "999"))
    ok &= run_modified("wrong t_star (0, in range)",
                       lambda b: find_cyc(b).__setitem__("t_star", 0))
    ok &= run_modified("wrong t_star (3, in range)",
                       lambda b: find_cyc(b).__setitem__("t_star", 3))
    ok &= run_modified("empty certificate set", lambda b: b.clear())
    ok &= run_modified(
        "missing K (drop one record)",
        lambda b: b.remove(next(r for r in b if r["K"] == 9)))
    ok &= run_modified(
        "duplicate K record",
        lambda b: b.append(dict(next(r for r in b if r["K"] == 9))))
    return ok


def main():
    ok = True
    ok &= t_params()
    ok &= t_scalar()
    ok &= t_word_detect()
    ok &= t_certificates()
    ok &= t_windows()
    ok &= t_alt4()
    ok &= t_lock()
    ok &= t_family_cover()
    ok &= t_theorem_b()
    ok &= t_theorem_c()
    ok &= t_theorem_d()
    ok &= t_d_independent_negative()
    ok &= t_bii_step_bound()
    ok &= t_wp4_error_semantics()
    ok &= t_a_plateau_allk()
    ok &= t_e_cycle_allk()
    ok &= t_b_lock_allk()
    ok &= t_c_alpha_endpoints()
    print(("ALL OK" if ok else "SOME FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
