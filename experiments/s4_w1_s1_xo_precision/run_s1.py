"""Runner for S4 Wave 1 / S1: XO high-precision replay.

One command executes the pre-registered round (config.json):

    python experiments/s4_w1_s1_xo_precision/run_s1.py

Stages (all gates frozen in config.json BEFORE any run):

  V  float64 regeneration: committed per-round CSVs are hash-checked
     against their recorded metadata, then the XO (Gate 2 runner code)
     and KC (Gate 1 runner code) trajectories are regenerated through
     the COMMITTED runner modules and every active round must be
     string-identical to the committed rows. Any mismatch halts the
     round before any precision work. The regenerated run objects (now
     proven identical to the committed record) provide full-precision
     float64 values for the numeric comparisons.
  R  precision ladder: exact Fraction tier (identity check on), mp100,
     mp50, XO at K in {50, 100, 200, 400}; KC sanity: exact @ K=50,
     mp100 @ all K. Feasibility gates per config; a gate failure skips
     the K for that tier and is recorded (never silently).
  C  comparisons against the F tier + itinerary analysis + metadata.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import io
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import mpmath as mp

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "experiments" / "s4_w1_s1_xo_precision"))

import numpy as np  # noqa: E402

import replay  # noqa: E402
from sand_m0.adaptive import run_m1a_deterministic  # noqa: E402
from sand_m0.model import (  # noqa: E402
    DISTRIBUTIONS,
    witness_opposed_distribution,
)

K_LIST = [50, 100, 200, 400]
OUT_DIR = Path(__file__).resolve().parent / "results"
GATE2_DIR = REPO_ROOT / "experiments" / "m2_gate2_sign_placement"
GATE1_DIR = REPO_ROOT / "experiments" / "m2_gate1_generality"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT,
            capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def csv_lines(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as fh:
        return fh.read().splitlines()


# --------------------------------------------------------------------------
# STAGE V — float64 regeneration and string verification
# --------------------------------------------------------------------------

def stage_v() -> tuple[dict, dict, list[dict]]:
    """Returns (gate2_integrity, kc_verification_info, verified float rows).

    float rows: list of dicts with kernel/K/t/j_t/d(float)/rebound(bool)/
    worst(float)/R(float)/boundary(bool)/d_after(float) taken from the
    regenerated (string-verified) run objects at full float precision.
    """
    print("=== STAGE V: float64 regeneration & string verification ===")
    v_rows = []

    # --- integrity of the committed reference CSVs ---
    integrity = {}
    for tag, exp_dir, meta_name, csv_name in [
        ("gate2", GATE2_DIR, "metadata.json", "per_round_diagnostics.csv"),
        ("gate1", GATE1_DIR, "metadata.json", "per_round_diagnostics.csv"),
    ]:
        meta = json.loads((exp_dir / "results" / meta_name).read_text(encoding="utf-8"))
        recorded = meta["outputs"][csv_name]
        actual = sha256_of(exp_dir / "results" / csv_name)
        integrity[tag] = {
            "recorded": recorded, "actual": actual, "match": recorded == actual}
        print(f"  [{tag}] committed CSV sha256 {'MATCHES' if recorded == actual else 'MISMATCHES'} metadata record")
        if recorded != actual:
            print("HALT: committed reference CSV does not match its recorded hash.")
            raise SystemExit(1)

    gate2 = load_module("gate2_runner", GATE2_DIR / "run_experiment.py")
    gate1 = load_module("gate1_runner", GATE1_DIR / "run_experiment.py")

    q_near = DISTRIBUTIONS["near"]
    xo = witness_opposed_distribution()
    alpha, tol, tau = 0.25, 1e-12, 1e-12
    float_rows: list[dict] = []

    # --- XO via the committed Gate-2 runner code path ---
    committed2 = list(csv.DictReader(
        (GATE2_DIR / "results" / "per_round_diagnostics.csv").open(encoding="utf-8")))
    fieldnames2 = csv_lines(
        GATE2_DIR / "results" / "per_round_diagnostics.csv")[0].split(",")
    for K in K_LIST:
        g = xo.bin_probs(K)
        run = run_m1a_deterministic(
            q_near, K=K, T_max=10 * K, alpha=alpha, stop_tolerance=tol,
            redistribution="throw", redistribution_dist=xo)
        per, _ = gate2.per_round_diagnostics("XO", K, run, g, xo, alpha,
                                             tau, 0.2)
        committed_k = [r for r in committed2
                       if r["kernel"] == "XO" and int(r["K"]) == K]
        mism = _compare_rows(per, committed_k, fieldnames2)
        n = len(per)
        v_rows.append({"source": f"XO K={K} (gate2 runner)", "rows": n,
                       "string_mismatches": mism,
                       "pass": mism == 0})
        print(f"  XO K={K}: {n} active rounds, {mism} string mismatches")
        if mism:
            print("HALT: regenerated XO rows are not string-identical.")
            raise SystemExit(1)
        states = np.asarray(run["states"], dtype=float)
        T = np.arange(1, K + 1) / K
        for t, row in enumerate(per):
            mass_t = states[t]
            D = np.cumsum(mass_t) - T
            float_rows.append({
                "kernel": "XO", "K": K, "t": t, "j_t": int(row["j_t"]),
                "d": float(row["D_max_before"]),
                "d_full": float(D[int(row["j_t"]) - 1]),
                "rebound": row["rebound"] is True,
                "worst": float(row["worst_excess"]),
                "R": float(row["R_t"]),
                "boundary": row["numerical_boundary"] is True,
                "moved": float(row["moved_mass"]),
            })
        d_after_final = float(run["rows"][n]["D_max"])
        float_rows[-1]["d_after"] = d_after_final

    # --- KC via the committed Gate-1 runner code path ---
    committed1 = list(csv.DictReader(
        (GATE1_DIR / "results" / "per_round_diagnostics.csv").open(encoding="utf-8")))
    fieldnames1 = csv_lines(
        GATE1_DIR / "results" / "per_round_diagnostics.csv")[0].split(",")
    kc_stop = {}
    for K in K_LIST:
        g = q_near.bin_probs(K)
        run = run_m1a_deterministic(
            q_near, K=K, T_max=10 * K, alpha=alpha, stop_tolerance=tol,
            redistribution="throw", redistribution_dist=None)
        per, _ = gate1.per_round_diagnostics("KC", K, run, g, alpha, tau, 0.2)
        committed_k = [r for r in committed1
                       if r["kernel"] == "KC" and int(r["K"]) == K]
        mism = _compare_rows(per, committed_k, fieldnames1)
        v_rows.append({"source": f"KC K={K} (gate1 runner)", "rows": len(per),
                       "string_mismatches": mism, "pass": mism == 0})
        print(f"  KC K={K}: {len(per)} active rounds, {mism} string mismatches")
        if mism:
            print("HALT: regenerated KC rows are not string-identical.")
            raise SystemExit(1)
        kc_stop[K] = {"t_stop": run["stopped_at"], "N_R": sum(
            1 for r in per if r["rebound"])}
        states = np.asarray(run["states"], dtype=float)
        T = np.arange(1, K + 1) / K
        for t, row in enumerate(per):
            mass_t = states[t]
            D = np.cumsum(mass_t) - T
            float_rows.append({
                "kernel": "KC", "K": K, "t": t, "j_t": int(row["j_t"]),
                "d": float(row["D_max_before"]),
                "d_full": float(D[int(row["j_t"]) - 1]),
                "rebound": row["rebound"] is True,
                "worst": float(row["worst_excess"]),
                "R": float(row["R_t"]),
                "boundary": row["numerical_boundary"] is True,
                "moved": float(row["moved_mass"]),
            })
        n = len(per)
        float_rows[-1]["d_after"] = float(run["rows"][n]["D_max"])

    write_csv(OUT_DIR / "stage_v_verification.csv",
              ["source", "rows", "string_mismatches", "pass"],
              [{**r, "pass": str(r["pass"])} for r in v_rows])
    print("  STAGE V PASSED — regenerated float64 rows are string-identical "
          "to the committed records.")
    return integrity, kc_stop, float_rows


def _compare_rows(regen: list[dict], committed: list[dict],
                  fieldnames: list[str]) -> int:
    """String-compare regenerated dicts against committed CSV rows through
    the same csv serialization (the committed writer's exact rendering)."""
    if len(regen) != len(committed):
        return -1
    mism = 0
    for rg, cm in zip(regen, committed):
        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fieldnames)
        writer.writerow({k: rg[k] for k in fieldnames})
        regen_line = buf.getvalue().rstrip("\r\n")
        comm_line = ",".join(cm[k] for k in fieldnames)
        if regen_line != comm_line:
            mism += 1
            if mism <= 2:
                print(f"    mismatch:\n      regen: {regen_line}\n      comm : {comm_line}")
    return mism


# --------------------------------------------------------------------------
# STAGE R — precision ladder
# --------------------------------------------------------------------------

def stage_r():
    print("\n=== STAGE R: precision ladder (exact / mp100 / mp50) ===")
    runs = {}
    tier_rows = []
    budgets = {50: 300, 100: 900, 200: 2700, 400: None}
    prev_ok = True
    for K in K_LIST:
        gate = budgets[K]
        if not prev_ok:
            print(f"  XO K={K} exact: SKIPPED (feasibility gate)")
            runs[("XO", K, "exact")] = None
            if gate is None:
                break
            continue
        r = replay.replay("XO", K, "exact", 10 * K, identity_check=True)
        halts = []
        if not r["exact_flags_ok"]:
            halts.append("exact invariant/selection-bound violation")
        if not r["identity_ok"]:
            halts.append("exact worst-identity violation")
        if halts:
            print(f"HALT for XO K={K}: {halts} — contradicts the registered "
                  "Model-derived result; interpretation paused.")
            raise SystemExit(1)
        within = gate is None or r["wall_seconds"] <= gate
        prev_ok = within
        print(f"  XO K={K} exact: {len(r['rows'])} rounds in {r['wall_seconds']}s "
              f"(gate {'pass' if within else 'FAIL'}), censored={r['censored']}, "
              f"exact_flags_ok={r['exact_flags_ok']}")
        if not within:
            runs[("XO", K, "exact")] = None
            continue
        runs[("XO", K, "exact")] = r
        for tier in ("mp100", "mp50"):
            runs[("XO", K, tier)] = replay.replay("XO", K, tier, 10 * K)
            print(f"  XO K={K} {tier}: {len(runs[('XO', K, tier)]['rows'])} rounds "
                  f"in {runs[('XO', K, tier)]['wall_seconds']}s")
    # KC sanity
    runs[("KC", 50, "exact")] = replay.replay("KC", 50, "exact", 500)
    print(f"  KC K=50 exact: t_stop={runs[('KC', 50, 'exact')]['stopped_at']} "
          f"in {runs[('KC', 50, 'exact')]['wall_seconds']}s")
    for K in K_LIST:
        runs[("KC", K, "mp100")] = replay.replay("KC", K, "mp100", 10 * K)
        print(f"  KC K={K} mp100: t_stop="
              f"{runs[('KC', K, 'mp100')]['stopped_at']} "
              f"in {runs[('KC', K, 'mp100')]['wall_seconds']}s")
    _write_tier_csvs(runs)
    return runs


FIELDS = ["t", "s", "d", "e", "moved", "worst", "R", "boundary",
          "tie_exact", "invariant_exact_ok", "s_bound_ok", "rebound"]


def _write_tier_csvs(runs: dict) -> None:
    for (kernel, K, tier), r in runs.items():
        if r is None:
            continue
        flags = replay.rebound_flags(r)
        rows = []
        for row, reb in zip(r["rows"], flags):
            rows.append({
                "t": row["t"], "s": row["s"], "d": row["d"], "e": row["e"],
                "moved": row["moved"], "worst": row["worst"], "R": row["R"],
                "boundary": str(bool(row["boundary"])),
                "tie_exact": (str(bool(row["tie_exact"]))
                              if row["tie_exact"] is not None else ""),
                "invariant_exact_ok": (str(bool(row["invariant_exact_ok"]))
                                       if row["invariant_exact_ok"] is not None else ""),
                "s_bound_ok": str(bool(row["s_bound_ok"])),
                "rebound": str(bool(reb)),
            })
        write_csv(OUT_DIR / f"{kernel.lower()}_{tier}_K{K}.csv", FIELDS, rows)


# --------------------------------------------------------------------------
# STAGE C — comparisons + itinerary + metadata
# --------------------------------------------------------------------------

def _tier_float(s: str) -> float:
    return float(s)


def stage_c(runs: dict, float_rows: list[dict], integrity: dict,
            kc_stop: dict) -> None:
    print("\n=== STAGE C: comparisons & itinerary analysis ===")
    summary_rows = []
    pair_rows = []
    itinerary_rows = []
    kc_rows = []

    for K in K_LIST:
        fx = [r for r in float_rows if r["kernel"] == "XO" and r["K"] == K]
        tiers = {}
        for tier in ("exact", "mp100", "mp50"):
            r = runs[("XO", K, tier)]
            if r is None:
                continue
            flags = replay.rebound_flags(r)
            tiers[tier] = {"rows": r["rows"], "flags": flags, "run": r}
        # --- pairwise comparisons vs F (and exact-vs-mp) ---
        f_s = [r["j_t"] for r in fx]
        f_reb = [r["rebound"] for r in fx]
        f_d_full = [r["d_full"] for r in fx]
        f_bnd = [r["boundary"] for r in fx]
        for tier, obj in tiers.items():
            t_s = [row["s"] for row in obj["rows"]]
            n = min(len(f_s), len(t_s))
            div = [i for i in range(n) if f_s[i] != t_s[i]]
            reb_t = obj["flags"][:n]
            reb_mism = [i for i in range(n) if f_reb[i] != reb_t[i]]
            d_tiers = [_tier_float(row["d"]) for row in obj["rows"]]
            max_dd = max(abs(a - b) for a, b in zip(f_d_full[:n], d_tiers[:n]))
            pair_rows.append({
                "K": K, "pair": f"F_vs_{tier}", "rounds_compared": n,
                "selection_first_divergence": div[0] if div else "",
                "selection_divergent_rounds": len(div),
                "rebound_flag_mismatches": len(reb_mism),
                "rebound_mismatch_rounds": ";".join(
                    str(i) for i in reb_mism[:20]) + ("..." if len(reb_mism) > 20 else ""),
                "max_abs_d_diff": f"{max_dd:.3e}",
            })
        if "mp100" in tiers and "mp50" in tiers:
            s100 = [r["s"] for r in tiers["mp100"]["rows"]]
            s50 = [r["s"] for r in tiers["mp50"]["rows"]]
            div = [i for i in range(min(len(s100), len(s50)))
                   if s100[i] != s50[i]]
            pair_rows.append({
                "K": K, "pair": "mp100_vs_mp50",
                "rounds_compared": min(len(s100), len(s50)),
                "selection_first_divergence": div[0] if div else "",
                "selection_divergent_rounds": len(div),
                "rebound_flag_mismatches": "", "rebound_mismatch_rounds": "",
                "max_abs_d_diff": "",
            })
        # --- per-tier summary ---
        n_f_active = len(fx)
        n_reb_f = sum(f_reb)
        for tag, obj in [("F_float64", {"rows": None, "flags": f_reb}),
                         ] + [(t, o) for t, o in tiers.items()]:
            if tag == "F_float64":
                s_vals = f_s
                flags = f_reb
                bnd = sum(f_bnd)
                ties = ""
                e_vals = None
                dvals = [r["d"] for r in fx]
                maxr = max(r["R"] for r in fx)
                wall = ""
                ex_ok = ""
            else:
                run = obj["run"]
                s_vals = [row["s"] for row in obj["rows"]]
                flags = obj["flags"]
                bnd = sum(1 for row in obj["rows"] if row["boundary"])
                ties = sum(1 for row in obj["rows"] if row["tie_exact"])
                e_vals = [_tier_float(row["e"]) for row in obj["rows"]]
                dvals = [_tier_float(row["d"]) for row in obj["rows"]]
                maxr = max(_tier_float(row["R"]) for row in obj["rows"])
                wall = run["wall_seconds"]
                ex_ok = str(run["exact_flags_ok"] and run["identity_ok"])
            n_reb = sum(flags)
            distinct = sorted(set(s_vals))
            m = (5 * K + 5) // 6
            row = {
                "kernel": "XO", "K": K, "tier": tag,
                "active_rounds": len(s_vals),
                "stop_status": ("exact_stop" if
                                (obj["run"] and obj["run"]["stopped_at"] is not None)
                                else "censored_at_H") if tag != "F_float64"
                else ("exact_stop" if False else "censored_at_H"),
                "N_R": n_reb, "f_R": f"{n_reb / max(len(flags), 1):.5f}",
                "boundary_rounds": bnd, "exact_ties": ties,
                "distinct_s": len(distinct),
                "s_max": max(distinct), "s_bound_m": m,
                "max_abs_e": f"{max(abs(v) for v in e_vals):.6e}" if e_vals else "",
                "final_d": f"{dvals[-1]:.12e}",
                "distinct_d_values": len(set(dvals)),
                "max_R_t": f"{maxr:.6f}",
                "exact_flags_ok": ex_ok, "wall_seconds": wall,
            }
            summary_rows.append(row)
        # --- itinerary (frozen scope: XO, all tiers) ---
        for tag, s_vals in ([("F_float64", f_s)]
                            + [(t, [row["s"] for row in o["rows"]])
                               for t, o in tiers.items()]):
            runs_list = []
            cur, ln = None, 0
            for s in s_vals:
                if s == cur:
                    ln += 1
                else:
                    if cur is not None:
                        runs_list.append((cur, ln))
                    cur, ln = s, 1
            runs_list.append((cur, ln))
            visits = {}
            for s in s_vals:
                visits[s] = visits.get(s, 0) + 1
            top_visits = sorted(visits.items(), key=lambda kv: -kv[1])[:10]
            trans = {}
            for a, b in zip(s_vals, s_vals[1:]):
                if a != b:
                    trans[(a, b)] = trans.get((a, b), 0) + 1
            top_trans = sorted(trans.items(), key=lambda kv: -kv[1])[:10]
            for s, visits_n in top_visits:
                rl = [ln for v, ln in runs_list if v == s]
                itinerary_rows.append({
                    "K": K, "tier": tag, "s": s, "visits": visits_n,
                    "runs": len(rl), "max_run": max(rl),
                    "mean_run": f"{sum(rl) / len(rl):.2f}",
                    "share": f"{visits_n / len(s_vals):.5f}",
                })
            it_last = itinerary_rows[-1]
            for (a, b), cnt in top_trans:
                itinerary_rows.append({
                    "K": K, "tier": tag, "s": f"{a}->{b}", "visits": cnt,
                    "runs": "", "max_run": "", "mean_run": "",
                    "share": f"{cnt / max(1, sum(trans.values())):.5f}",
                })

    # --- KC comparison ---
    for K in K_LIST:
        row = {"K": K, "t_stop_F": kc_stop[K]["t_stop"],
               "N_R_F": kc_stop[K]["N_R"]}
        r100 = runs[("KC", K, "mp100")]
        row["t_stop_mp100"] = r100["stopped_at"]
        row["N_R_mp100"] = sum(replay.rebound_flags(r100))
        if K == 50:
            rx = runs[("KC", 50, "exact")]
            row["t_stop_exact"] = rx["stopped_at"]
            row["N_R_exact"] = sum(replay.rebound_flags(rx))
        else:
            row["t_stop_exact"] = ""
            row["N_R_exact"] = ""
        kc_rows.append(row)
        print(f"  KC K={K}: t_stop F={row['t_stop_F']} "
              f"mp100={row['t_stop_mp100']} exact={row['t_stop_exact']}")

    write_csv(OUT_DIR / "comparison_pairs.csv", list(pair_rows[0]), pair_rows)
    write_csv(OUT_DIR / "summary_by_tier.csv", list(summary_rows[0]), summary_rows)
    write_csv(OUT_DIR / "itinerary_xo.csv",
              ["K", "tier", "s", "visits", "runs", "max_run", "mean_run", "share"],
              itinerary_rows)
    write_csv(OUT_DIR / "comparison_kc.csv", list(kc_rows[0]), kc_rows)

    print("\n--- summary (XO, per K / tier) ---")
    for r in summary_rows:
        print(f"  K={r['K']:>3} {r['tier']:<9} N_R={r['N_R']:>4} "
              f"f_R={r['f_R']} bnd/tie={r['boundary_rounds']}/{r['exact_ties']} "
              f"distinct_s={r['distinct_s']} s_max={r['s_max']} "
              f"|e|max={r['max_abs_e'] or '-'} dvals={r['distinct_d_values']}")
    print("\n--- pairwise (selection / rebound / d) ---")
    for r in pair_rows:
        print(f"  K={r['K']:>3} {r['pair']:<14} sel_div={r['selection_divergent_rounds']} "
              f"first={r['selection_first_divergence']} "
              f"reb_mism={r['rebound_flag_mismatches']} "
              f"max|dd|={r['max_abs_d_diff']}")


def main() -> int:
    t0 = time.time()
    cfg = json.loads((Path(__file__).resolve().parent / "config.json")
                     .read_text(encoding="utf-8"))
    OUT_DIR.mkdir(exist_ok=True)
    integrity, kc_stop, float_rows = stage_v()
    runs = stage_r()
    stage_c(runs, float_rows, integrity, kc_stop)
    output_files = sorted(p.name for p in OUT_DIR.iterdir()
                          if p.is_file() and p.name != "metadata.json")
    metadata = {
        "model": cfg["model"], "model_version": cfg["model_version"],
        "experiment": cfg["experiment"],
        "run_date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": cfg,
        "gate2_integrity": integrity,
        "environment": {"python": sys.version, "platform": platform.platform(),
                        "numpy": np.__version__, "mpmath": mp.__version__},
        "git": {"commit": git_commit()},
        "runtime_seconds": round(time.time() - t0, 2),
        "command": "python experiments/s4_w1_s1_xo_precision/run_s1.py",
        "outputs": {name: sha256_of(OUT_DIR / name) for name in output_files},
    }
    (OUT_DIR / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\noutputs written to {OUT_DIR} "
          f"(runtime {round(time.time() - t0, 2)} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
