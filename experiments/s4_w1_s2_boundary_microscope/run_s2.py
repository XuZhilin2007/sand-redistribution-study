"""Runner for S4 Wave 1 / S2: XO Boundary Switching Microscope.

Executes the frozen experiment (config.json, committed e706a0d — never
modified in-round):

    python experiments/s4_w1_s2_boundary_microscope/run_s2.py

  Step 0  committed-CSV gate (S0a-S0e): analysis-only pre-check on the
          committed S1 exact CSVs. CSV-derived values carry 60-digit
          storage truncation -> tolerance 1e-50 (frozen); failures halt
          before the census.
  Step 1  exact K census, K = 6..53 + anchors {100, 200, 400}, H = 10K,
          exact Fraction arithmetic (primary and only evidence tier),
          frozen classification (t* / R_m >= K / CYCLE / OTHER subtypes),
          mod-6 prediction comparison.
  Step 2  reduction verification C1-C7 per K (in-memory checks are EXACT
          Fraction equalities).

Hard halts (H1-H4) stop everything; scientific events (E1-E4) are
recorded and reported, never repaired in-round; no in-round scope
expansion is possible.
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

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "s4_w1_s1_xo_precision"))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "s4_w1_s2_boundary_microscope"))

import mpmath as mp  # noqa: E402

import replay  # noqa: E402  (S1 exact engine, six-layer verified)

ALPHA = Fraction(1, 4)
TOL = Fraction(1, 10 ** 12)
CENSUS = list(range(6, 54))
ANCHORS = [100, 200, 400]
ALL_K = CENSUS + ANCHORS
STEP0_TOL = Fraction(1, 10 ** 50)

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
S1_RESULTS = REPO_ROOT / "experiments" / "s4_w1_s1_xo_precision" / "results"
S1_CSV_COLUMNS = ["t", "s", "d", "e", "moved", "worst", "R", "boundary",
                  "tie_exact", "invariant_exact_ok", "s_bound_ok", "rebound"]

HARD = []          # hard-halt reasons
EVENTS = []        # scientific events E1-E4


def hard_halt(reason: str):
    HARD.append(reason)
    raise HardHalt(reason)


class HardHalt(Exception):
    pass


def event(kind: str, K: int, detail: str):
    EVENTS.append({"kind": kind, "K": K, "detail": detail})


def frac_of(dec: str) -> Fraction:
    """Exact Fraction from a decimal string (exact for terminating decimals)."""
    return Fraction(dec)


def sha256_of(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            d.update(chunk)
    return d.hexdigest()


def git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


# --------------------------------------------------------------------------
# frozen classification (pure function; unit-tested)
# --------------------------------------------------------------------------

def classify_itinerary(s_list: list[int], m: int, K: int) -> dict:
    """Frozen classifier (config.classification_frozen).

    t*   = first entry into {m-1, m} (None if never)
    R_m  = length of the final consecutive run of s = m over the full
           itinerary (0 if s does not end at m)
    LOCK      iff R_m >= K
    CYCLE     else iff post-entry set subset of {m-1, m} and both occur
    OTHER     else, subtype among: never_enters | interior_return |
              m_minus_1_only | short_tail_lock | other
    """
    t_star = next((t for t, s in enumerate(s_list) if s in (m - 1, m)), None)
    R_m = 0
    for s in reversed(s_list):
        if s == m:
            R_m += 1
        else:
            break
    if t_star is None:
        return {"t_star": None, "R_m": R_m, "class": "OTHER",
                "subtype": "never_enters", "post_set": []}
    post = s_list[t_star:]
    post_set = sorted(set(post))
    if R_m >= K:
        cls, subtype = "LOCK", ""
    elif set(post) <= {m - 1, m} and len(post_set) == 2:
        cls, subtype = "CYCLE", ""
    elif set(post) == {m}:
        cls, subtype = "OTHER", "short_tail_lock"
    elif set(post) == {m - 1}:
        cls, subtype = "OTHER", "m_minus_1_only"
    elif min(post) < m - 1:
        cls, subtype = "OTHER", "interior_return"
    else:
        cls, subtype = "OTHER", "other"
    return {"t_star": t_star, "R_m": R_m, "class": cls, "subtype": subtype,
            "post_set": post_set}


def tie_tail_length(raw_rows: list[dict]) -> int:
    """Final consecutive run of exact-tie rounds (secondary criterion)."""
    n = 0
    for r in reversed(raw_rows):
        if r["tie_exact"]:
            n += 1
        else:
            break
    return n


# --------------------------------------------------------------------------
# STEP 0 — committed-CSV gate (S0a-S0e)
# --------------------------------------------------------------------------

def mod6_predicted_class(K: int) -> str:
    return "LOCK" if K % 6 in (0, 1, 2) else "CYCLE"


def step0() -> list[dict]:
    print("=== STEP 0: committed-CSV gate (S0a-S0e, tolerance 1e-50) ===")
    out_rows = []
    all_ok = True
    for K in (50, 100, 200, 400):
        m, b = replay.anchor_m_b(K)
        Fm = 2 * b - b * b
        qK = Fraction(6 - (K % 6), 3 * K)
        rows = list(csv.DictReader(
            (S1_RESULTS / f"xo_exact_K{K}.csv").open(encoding="utf-8")))
        s = [int(r["s"]) for r in rows]
        # S0a classification
        cls = classify_itinerary(s, m, K)
        s1_registered = {50: "LOCK", 100: "CYCLE", 200: "LOCK", 400: "CYCLE"}[K]
        ok_a = cls["class"] == s1_registered
        # S0b F_m/4 == moved on s=m rounds
        bad_b = [t for t, r in enumerate(rows) if s[t] == m
                 and abs(frac_of(r["moved"]) - Fm / 4) > STEP0_TOL]
        ok_b = not bad_b
        # S0c lock signature e == 0 on s=m rounds, t >= 3 (K=50/200)
        if K in (50, 200):
            bad_c = [t for t, r in enumerate(rows) if s[t] == m and t >= 3
                     and frac_of(r["e"]) != 0]
        else:
            bad_c = []
        ok_c = not bad_c
        # S0d decomposition identities (exact integers)
        ties = sum(1 for t in range(cls["t_star"], len(rows))
                   if rows[t]["tie_exact"] == "True")
        rebs = sum(1 for t in range(cls["t_star"], len(rows))
                   if rows[t]["rebound"] == "True")
        cons = sum(1 for r in rows
                   if r["rebound"] == "False" and r["tie_exact"] == "False")
        A = sum(1 for x in s[cls["t_star"]:] if x == m)
        B = sum(1 for x in s[cls["t_star"]:] if x == m - 1)
        ok_d = (A == ties + rebs) and (cons == B + cls["t_star"])
        # S0e branch-B recovery
        p_rec = {t: Fraction(1, K) - frac_of(r["e"])
                 for t, r in enumerate(rows)
                 if s[t] == m - 1 and frac_of(r["e"]) > 0}
        bb_pairs = bb_ok = ba_cross = ba_ok = 0
        worst_resid = Fraction(0)
        for t, p_t in p_rec.items():
            if t + 1 in p_rec:
                bb_pairs += 1
                resid = abs((p_rec[t + 1] - p_t)
                            - Fraction(1, 4) * (Fm - p_t) * qK)
                worst_resid = max(worst_resid, resid)
                bb_ok += int(resid <= STEP0_TOL)
            elif t + 1 < len(rows) and s[t + 1] == m:
                ba_cross += 1
                ba_ok += int(p_t + Fraction(1, 4) * (Fm - p_t) * qK
                             > Fraction(1, K))
        ok_e = (bb_ok == bb_pairs) and (ba_ok == ba_cross)
        row_ok = ok_a and ok_b and ok_c and ok_d and ok_e
        all_ok = all_ok and row_ok
        out_rows.append({
            "K": K, "S0a_class": cls["class"], "S0a_ok": ok_a,
            "S0b_bad_rounds": len(bad_b), "S0b_ok": ok_b,
            "S0c_bad_rounds": len(bad_c), "S0c_ok": ok_c,
            "S0d_A_vs_ties_rebs": f"{A}=={ties}+{rebs}",
            "S0d_cons_vs_B_tstar": f"{cons}=={B}+{cls['t_star']}",
            "S0d_ok": ok_d,
            "S0e_bb_pairs": f"{bb_ok}/{bb_pairs}", "S0e_ba_cross": f"{ba_ok}/{ba_cross}",
            "S0e_worst_residual": f"{float(worst_resid):.3e}", "S0e_ok": ok_e,
            "pass": row_ok,
        })
        print(f"  K={K}: class={cls['class']} (S0a {ok_a}) S0b {ok_b} "
              f"S0c {ok_c} S0d {ok_d} S0e bb={bb_ok}/{bb_pairs} "
              f"ba={ba_ok}/{ba_cross} resid={float(worst_resid):.2e} -> "
              f"{'PASS' if row_ok else 'FAIL'}")
    write_csv(OUT / "step0_s1_csv_checks.csv", list(out_rows[0]), out_rows)
    if not all_ok:
        hard_halt("Step-0 gate failed (S0a/S0b/S0c are hard halts; "
                  "S0d/S0e failures halt pending review)")
    print("  STEP 0 PASSED\n")
    return out_rows


# --------------------------------------------------------------------------
# STEP 1 + 2 — exact census with reduction verification
# --------------------------------------------------------------------------

REDUCTION_FIELDS = ["t", "s", "branch", "p_m", "F_m", "restriction_ok",
                    "branch_eq_ok", "sel_branch_ok", "dpred_ok",
                    "class_predicted", "class_observed", "tie_exact"]


def exact_replay_s2(K: int, instrument_full: bool, H: int | None = None) -> dict:
    """Extended exact replay: S1 semantics + boundary diagnostics.

    In-memory check values are EXACT Fractions. Shared-column S1 rows are
    produced only for instrument_full (anchors, H2 string check).
    H defaults to the frozen horizon 10*K; tests use short slices.
    """
    t0 = time.perf_counter()
    H = 10 * K if H is None else H
    m, b = replay.anchor_m_b(K)
    anchor = b * (1 - b)
    Fm_const = 2 * b - b * b
    one_over_K = Fraction(1, K)
    zero = Fraction(0)
    p = replay.exact_bin_probs(replay.near_cdf_exact, K)
    q = replay.exact_bin_probs(replay.xo_cdf_exact, K)
    qK = q[m - 1]
    raw = []
    s_list = []
    d_list = []
    c1_bad_states = []
    # C1 state 0
    if sum(p[:m]) != Fm_const:
        c1_bad_states.append(0)
    for t in range(H):
        D = replay.cumsum_excess(p, K)
        s = replay.argmax_first(D)
        d = D[s - 1]
        if d <= TOL:
            hard_halt(f"K={K}: exact stop observed at t={t} — contradicts "
                      "the non-stopping consequence of [Theorem]")
        p_m = p[m - 1]
        F_m = sum(p[:m])
        interior_max = max(D[0:m - 2])
        pair_max = max(D[m - 2], D[m - 1])
        restriction_ok = interior_max < pair_max
        p_m_tie = (D[m - 2] == D[m - 1])
        branch = "A" if s == m else ("B" if s == m - 1 else "pre")
        sel_ok = ((p_m > one_over_K) if s == m else (p_m <= one_over_K))
        dpred = anchor + (one_over_K - p_m if p_m < one_over_K else zero)
        dpred_ok = (d == dpred)
        if branch == "A":
            p_pred = Fraction(3, 4) * p_m + Fraction(1, 4) * F_m * qK
        elif branch == "B":
            p_pred = p_m + Fraction(1, 4) * (F_m - p_m) * qK
        else:
            p_pred = None
        if instrument_full:
            dg = replay.diagnostics(p, q, K, s, ALPHA)
            worst, R, moved = dg["worst"], dg["R"], dg["moved"]
        else:
            worst = R = moved = None
        raw.append({
            "t": t, "s": s, "branch": branch,
            "p_m": p_m, "F_m": F_m, "_d": d, "_worst": worst, "_R": R,
            "_moved": moved,
            "restriction_ok": restriction_ok, "_p_m_tie": p_m_tie,
            "_p_pred": p_pred, "sel_branch_ok": sel_ok, "dpred_ok": dpred_ok,
            "invariant_ok": (D[m - 1] == anchor), "s_bound_ok": (s <= m),
        })
        s_list.append(s)
        d_list.append(d)
        p = replay.update_state(p, q, s, ALPHA)
        if sum(p[:m]) != Fm_const:
            c1_bad_states.append(t + 1)
    n = len(raw)
    # post-hoc: d increments (incl. final state), branch equations, classes
    D_final = replay.cumsum_excess(p, K)
    d_final = D_final[replay.argmax_first(D_final) - 1]
    p_m_final = p[m - 1]
    for i, r in enumerate(raw):
        d_next = d_list[i + 1] if i + 1 < n else d_final
        r["_d_next"] = d_next
        r["tie_exact"] = (d_next == r["_d"])
        r["class_observed"] = ("tie" if d_next == r["_d"]
                               else ("rebound" if d_next > r["_d"]
                                     else "contraction"))
        if r["_p_pred"] is not None:
            p_m_next = (raw[i + 1]["p_m"] if i + 1 < n else p_m_final)
            r["branch_eq_ok"] = (p_m_next == r["_p_pred"])
        else:
            r["branch_eq_ok"] = None
        # frozen round-class prediction
        if r["branch"] == "A":
            r["class_predicted"] = ("tie" if p_m_next >= one_over_K
                                    else "rebound")
        elif r["branch"] == "B":
            r["class_predicted"] = "contraction"
        else:
            r["class_predicted"] = ""
    cls = classify_itinerary(s_list, m, K)
    # C6 reduced orbit (applicable iff C2 & C3 hold post-entry)
    first_restr = next((r["t"] for r in raw[cls["t_star"]:]
                        if not r["restriction_ok"]), None) if cls["t_star"] is not None else None
    first_brancheq = next((r["t"] for r in raw[cls["t_star"]:]
                           if r["branch_eq_ok"] is False), None) if cls["t_star"] is not None else None
    c6_applicable = (cls["t_star"] is not None and first_restr is None
                     and first_brancheq is None)
    c6_ok = None
    first_c6 = None
    if c6_applicable:
        p_map = raw[cls["t_star"]]["p_m"]
        c6_ok = True
        for i in range(cls["t_star"], n):
            r = raw[i]
            if p_map > one_over_K:
                p_map = Fraction(3, 4) * p_map + Fraction(1, 4) * Fm_const * qK
                pred = ("tie" if p_map >= one_over_K else "rebound")
            else:
                p_map = p_map + Fraction(1, 4) * (Fm_const - p_map) * qK
                pred = "contraction"
            p_full_next = (raw[i + 1]["p_m"] if i + 1 < n else p_m_final)
            if p_map != p_full_next or pred != r["class_observed"]:
                c6_ok = False
                first_c6 = i
                break
    return {
        "K": K, "m": m, "b": b, "qK": qK, "H": H, "raw": raw, "cls": cls,
        "c1_bad_states": c1_bad_states, "c6_applicable": c6_applicable,
        "c6_ok": c6_ok, "first_c6": first_c6,
        "first_restr_violation": first_restr,
        "first_branch_eq_mismatch": first_brancheq,
        "n_post_restriction_violations":
            sum(1 for r in raw[cls["t_star"]:] if not r["restriction_ok"])
            if cls["t_star"] is not None else 0,
        "n_p_m_tie_rounds": sum(1 for r in raw if r["_p_m_tie"]),
        "first_sel_mismatch": next((r["t"] for r in raw[cls["t_star"]:]
                                    if not r["sel_branch_ok"]),
                                   None) if cls["t_star"] is not None else None,
        "first_dpred_mismatch": next((r["t"] for r in raw[cls["t_star"]:]
                                      if not r["dpred_ok"]),
                                     None) if cls["t_star"] is not None else None,
        "tie_tail": tie_tail_length(raw),
        "N_R": sum(1 for r in raw if r["class_observed"] == "rebound"),
        "wall_seconds": round(time.perf_counter() - t0, 2),
        "final_p_m": p_m_final,
    }


def s1_shared_row(res: dict, i: int) -> dict:
    """S1-format shared columns (anchor H2 check)."""
    r = res["raw"][i]
    anchor = res["b"] * (1 - res["b"])
    return {
        "t": r["t"], "s": r["s"], "d": replay._to_storage(r["_d"]),
        "e": replay._to_storage(r["_d"] - anchor),
        "moved": replay._to_storage(r["_moved"]),
        "worst": replay._to_storage(r["_worst"]),
        "R": replay._to_storage(r["_R"]),
        "boundary": str(bool(abs(r["_worst"]) <= TOL)),
        "tie_exact": str(bool(r["tie_exact"])),
        "invariant_exact_ok": str(bool(r["invariant_ok"])),
        "s_bound_ok": str(bool(r["s_bound_ok"])),
        "rebound": str(r["class_observed"] == "rebound"),
    }


def h2_anchor_check(res: dict) -> None:
    """H2: shared columns must string-match the committed S1 exact CSV."""
    K = res["K"]
    committed = list(csv.DictReader(
        (S1_RESULTS / f"xo_exact_K{K}.csv").open(encoding="utf-8")))
    if len(committed) != len(res["raw"]):
        hard_halt(f"H2 K={K}: row count {len(res['raw'])} != committed {len(committed)}")
    for i, cr in enumerate(committed):
        mine = s1_shared_row(res, i)
        for col in S1_CSV_COLUMNS:
            if str(mine[col]) != cr[col]:
                hard_halt(f"H2 K={K}: shared-column mismatch at t={i}, "
                          f"column {col}: '{mine[col]}' != '{cr[col]}'")
    reg = {50: "LOCK", 100: "CYCLE", 200: "LOCK", 400: "CYCLE"}[K]
    if res["cls"]["class"] != reg:
        hard_halt(f"H2 K={K}: class {res['cls']['class']} != S1-registered {reg}")


def reduction_rows(res: dict) -> list[dict]:
    rows = []
    for r in res["raw"]:
        rows.append({
            "t": r["t"], "s": r["s"], "branch": r["branch"],
            "p_m": replay._to_storage(r["p_m"]),
            "F_m": replay._to_storage(r["F_m"]),
            "restriction_ok": str(bool(r["restriction_ok"])),
            "branch_eq_ok": ("" if r["branch_eq_ok"] is None
                             else str(bool(r["branch_eq_ok"]))),
            "sel_branch_ok": str(bool(r["sel_branch_ok"])),
            "dpred_ok": str(bool(r["dpred_ok"])),
            "class_predicted": r["class_predicted"],
            "class_observed": r["class_observed"],
            "tie_exact": str(bool(r["tie_exact"])),
        })
    return rows


def main() -> int:
    t_start = time.time()
    OUT.mkdir(exist_ok=True)
    cfg_hash = sha256_of(Path(__file__).resolve().parent / "config.json")
    print(f"S2 Boundary Switching Microscope — frozen config sha256 {cfg_hash[:16]}…")
    print(f"HEAD {git_commit()}; exact tier only; H = 10K; "
          f"census 6..53 + anchors {ANCHORS}\n")

    step0()

    print(f"=== STEP 1+2: exact census & reduction verification ===")
    census_rows = []
    for K in ALL_K:
        instrument_full = K in ANCHORS
        res = exact_replay_s2(K, instrument_full)
        if res["c1_bad_states"]:
            hard_halt(f"H3 K={K}: C1 invariant violated at states "
                      f"{res['c1_bad_states'][:5]} (contradicts [Theorem])")
        if K in ANCHORS:
            h2_anchor_check(res)
        predicted = mod6_predicted_class(K)
        match = (res["cls"]["class"] == predicted)
        if not match:
            event("E3", K, f"class {res['cls']['class']} != predicted {predicted}")
        if res["cls"]["class"] == "OTHER":
            event("E4", K, f"subtype {res['cls']['subtype']}")
        if res["first_branch_eq_mismatch"] is not None:
            event("E1", K, f"branch equation mismatch first at t="
                  f"{res['first_branch_eq_mismatch']}")
        if res["first_restr_violation"] is not None:
            event("E2", K, f"restriction violation first at t="
                  f"{res['first_restr_violation']}")
        write_csv(OUT / f"reduction_K{K}.csv", REDUCTION_FIELDS,
                  reduction_rows(res))
        desc = ">".join(str(r["s"]) for r in
                        res["raw"][:res["cls"]["t_star"]] if True) or "-"
        census_rows.append({
            "K": K, "m": res["m"], "b": replay._to_storage(res["b"]),
            "qK": f"{res['qK'].numerator}/{res['qK'].denominator}",
            "qK_closed_form_ok": str(res["qK"] == Fraction(6 - (K % 6), 3 * K)),
            "t_star": res["cls"]["t_star"],
            "descent_itinerary": desc,
            "R_m": res["cls"]["R_m"], "tie_tail": res["tie_tail"],
            "class_observed": res["cls"]["class"],
            "subtype": res["cls"]["subtype"],
            "class_predicted_mod6": predicted,
            "prediction_match": match,
            "N_R_exact": res["N_R"],
            "C1_ok": not res["c1_bad_states"],
            "C2_first_violation": res["first_restr_violation"] if res["first_restr_violation"] is not None else "",
            "C3_first_mismatch": res["first_branch_eq_mismatch"] if res["first_branch_eq_mismatch"] is not None else "",
            "C4_first_mismatch": res["first_sel_mismatch"] if res["first_sel_mismatch"] is not None else "",
            "C5_first_mismatch": res["first_dpred_mismatch"] if res["first_dpred_mismatch"] is not None else "",
            "C6_applicable": res["c6_applicable"],
            "C6_ok": "" if res["c6_ok"] is None else res["c6_ok"],
            "C6_first_mismatch": res["first_c6"] if res["first_c6"] is not None else "",
            "n_post_restriction_violations": res["n_post_restriction_violations"],
            "n_p_m_tie_rounds": res["n_p_m_tie_rounds"],
            "wall_seconds": res["wall_seconds"],
        })
        print(f"  K={K:>3}: {res['cls']['class']:<5} (pred {predicted:<5} "
              f"{'MATCH' if match else 'MISMATCH'}) t*={res['cls']['t_star']} "
              f"R_m={res['cls']['R_m']} N_R={res['N_R']} C6="
              f"{('ok' if res['c6_ok'] else 'FAIL') if res['c6_applicable'] else 'n/a'} "
              f"[{res['wall_seconds']}s]")
    write_csv(OUT / "per_k_census.csv", list(census_rows[0]), census_rows)

    n_match = sum(1 for r in census_rows if r["prediction_match"])
    n_c6 = sum(1 for r in census_rows if r["C6_applicable"])
    n_c6_ok = sum(1 for r in census_rows if r["C6_ok"] is True)
    print(f"\n=== SUMMARY ===")
    print(f"mod-6 prediction: {n_match}/{len(census_rows)} MATCH")
    print(f"reduction C6: applicable {n_c6}, exact-orbit ok {n_c6_ok}")
    print(f"scientific events: {len(EVENTS)}")
    for ev in EVENTS[:20]:
        print(f"  {ev['kind']} K={ev['K']}: {ev['detail']}")

    output_files = sorted(p.name for p in OUT.iterdir()
                          if p.is_file() and p.name != "metadata.json")
    metadata = {
        "model": "S4 Wave 1 / S2 — XO Boundary Switching Microscope",
        "run_date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "frozen_config_sha256": cfg_hash,
        "git": {"commit": git_commit()},
        "halted": False,
        "scientific_events": EVENTS,
        "summary": {"prediction_match": f"{n_match}/{len(census_rows)}",
                    "C6_applicable": n_c6, "C6_ok": n_c6_ok},
        "environment": {"python": sys.version,
                        "mpmath": mp.__version__},
        "runtime_seconds": round(time.time() - t_start, 2),
        "command": "python experiments/s4_w1_s2_boundary_microscope/run_s2.py",
        "outputs": {name: sha256_of(OUT / name) for name in output_files},
    }
    (OUT / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\noutputs in {OUT} (total {round(time.time() - t_start, 1)} s)")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except HardHalt as exc:
        print(f"\nHARD HALT: {exc}")
        (OUT / "HALT.txt").write_text(
            f"HARD HALT at {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n"
            f"HEAD {git_commit()}\n{exc}\n", encoding="utf-8")
        raise SystemExit(1)
