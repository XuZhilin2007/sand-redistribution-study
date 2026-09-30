"""Runner for M1C.1: fixed-uniform kernel ablation (Variant A vs U vs B).

One command regenerates every CSV, figure and the metadata file:

    python experiments/m1c1_uniform_kernel_ablation/run_experiment.py

Research question (frozen): if the M1 selection rule remains unchanged and
redistribution remains fixed and deficit-unaware, does removing the
near-bias eliminate or substantially weaken churn?

  * Variant A — canonical M1 baseline, untouched: initial state q_near,
    removed mass resprayed with fixed q_near.
  * Variant U — the M1C.1 middle control: same selection, same initial
    state (q_near), removed mass resprayed with the fixed uniform law
    q_uniform[i] = 1/K. Still fixed, state-unaware, deficit-unaware; only
    the near-bias is removed. Walks the same canonical fixed-throw code
    path as A (redistribution="throw", redistribution_dist=uniform).
  * Variant B — existing M1C oracle control: target-aware deficit filling.

A and B are re-run under the frozen configuration and string-compared
against the committed M1C summary (baseline integrity gate) before any
A/U/B interpretation. Variant U is the only new mechanism; no post-hoc
tuning of any kind.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from sand_m0.adaptive import run_m1a_deterministic  # noqa: E402
from sand_m0.model import DISTRIBUTIONS  # noqa: E402
from sand_m0.plotting import (  # noqa: E402
    plot_aub_trajectories,
    plot_aub_kernel_windows,
)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_state() -> dict:
    def run(args: list[str]) -> str:
        try:
            return subprocess.run(
                ["git", *args], cwd=REPO_ROOT, capture_output=True, text=True,
                check=True,
            ).stdout.strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            return "unknown"

    return {
        "commit": run(["rev-parse", "HEAD"]),
        "working_tree_dirty_at_run_time": bool(run(["status", "--porcelain"])),
    }


def read_csv_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def band(v: np.ndarray) -> dict:
    return {"p5": float(np.percentile(v, 5)), "p50": float(np.percentile(v, 50)),
            "p95": float(np.percentile(v, 95))}


def run_variant(variant: str, K: int, *, q_near, q_unif, T_max, alpha, tol) -> dict:
    """One frozen-configuration run; U reuses A's fixed-throw code path."""
    kwargs = {"redistribution": "throw"}
    if variant == "U":
        kwargs["redistribution_dist"] = q_unif
    elif variant == "B":
        kwargs = {"redistribution": "deficit_fill"}
    run = run_m1a_deterministic(q_near, K=K, T_max=T_max, alpha=alpha,
                                stop_tolerance=tol, **kwargs)
    stop = run["stopped_at"]
    n = stop if stop is not None else T_max
    rows = run["rows"]
    states = np.asarray(run["states"], dtype=float)
    Ud = np.array([K * r["U_L2"] for r in rows])
    Dm = np.array([r["D_max"] for r in rows[:n]])
    moved = np.array([r["removed_mass"] for r in rows[:n]])
    respray = q_unif if variant == "U" else q_near
    F_r = np.array([float(respray.cdf(np.array([r["a_t"]]))[0]) for r in rows[:n]])
    L1 = np.abs(states[1:n + 1] - states[:n]).sum(axis=1)
    return {"run": run, "stop": stop, "n": n, "Ud": Ud, "Dm": Dm, "moved": moved,
            "L1": L1, "F_r": F_r, "states": states}


def summary_row(variant: str, K: int, v: dict, alpha: float) -> dict:
    stop, Ud, Dm, moved = v["stop"], v["Ud"], v["Dm"], v["moved"]
    n = v["n"]
    traj_Ud = Ud[:n + 1]
    u0, u_best = float(traj_Ud[0]), float(traj_Ud.min())
    t_best = int(traj_Ud.argmin())
    u_stop = float(traj_Ud[-1])
    Db, Ub = band(Dm), band(traj_Ud)
    last25 = slice(max(int(0.75 * n), 0), n)
    late_band = Dm < 0.005
    dU = np.diff(traj_Ud)
    row = {
        "variant": variant, "K": K,
        "t_exact_stop": stop if stop is not None else "none",
        "active_rounds": n,
        "t_U_best": t_best,
        "U_density_initial": f"{u0:.6e}",
        "U_density_best": f"{u_best:.6e}",
        "U_density_at_stop": f"{u_stop:.6e}",
        "improvement_retention": f"{1.0 - (u_stop - u_best) / (u0 - u_best):.6f}",
        "D_max_p5": f"{Db['p5']:.5f}", "D_max_p50": f"{Db['p50']:.5f}",
        "D_max_p95": f"{Db['p95']:.5f}",
        "U_density_p5": f"{Ub['p5']:.5f}", "U_density_p50": f"{Ub['p50']:.5f}",
        "U_density_p95": f"{Ub['p95']:.5f}",
        "moved_mass_median": f"{float(np.median(moved)):.6f}",
        "moved_mass_median_last25pct": f"{float(np.median(moved[last25])):.6f}",
        "moved_mass_median_late_band":
            f"{float(np.median(moved[late_band])):.6f}" if late_band.any() else "n/a",
        "moved_mass_total": f"{float(moved.sum()):.6f}",
        "share_delta_U_positive": f"{float((dU > 0).mean()):.4f}",
        "share_D_max_increases": f"{float((np.diff(Dm) > 0).mean()):.4f}",
        "near_mass_05_at_stop": f"{v['run']['rows'][n]['near_mass_05']:.6f}",
    }
    return row


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=Path(__file__).resolve().parent / "config.json",
    )
    args = parser.parse_args(argv)

    config_path = args.config.resolve()
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    out_dir = config_path.parent / cfg["output_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)

    K_list = list(cfg["K_list"])
    T_max = cfg["T_max"]
    alpha = cfg["alpha"]
    tol = cfg["stop_tolerance_exact"]
    eps_list = list(cfg["epsilon_overlay"])
    q_near, q_unif = DISTRIBUTIONS["near"], DISTRIBUTIONS["uniform"]
    started = time.time()

    print(f"M1C.1 fixed-uniform kernel ablation — config: {config_path}")
    print(f"  K={K_list}, alpha={alpha}, exact tol={tol}, T_max={T_max}")

    runs: dict[tuple[str, int], dict] = {}
    for variant in ["A", "U", "B"]:
        for K in K_list:
            runs[(variant, K)] = run_variant(variant, K, q_near=q_near,
                                             q_unif=q_unif, T_max=T_max,
                                             alpha=alpha, tol=tol)
        print(f"  [{variant}] runs done: stops = "
              f"{[runs[(variant, K)]['stop'] for K in K_list]}")

    # summary rows for all variants (A/B rows double as the integrity cache)
    summary_rows = [summary_row(variant, K, runs[(variant, K)], alpha)
                    for variant in ["A", "U", "B"] for K in K_list]
    cache = {(r["variant"], int(r["K"])): r for r in summary_rows}

    # U invariant verification (per round, from states)
    for K in K_list:
        v = runs[("U", K)]
        states, rows, stop = v["states"], v["run"]["rows"], v["stop"]
        worst = {"mass_err": 0.0, "min_bin": np.inf,
                 "add_const_err": 0.0, "redep_err": 0.0}
        for t in range(stop):
            r = rows[t]
            p_minus = states[t].copy()
            p_minus[: r["j_t"]] *= (1.0 - alpha)
            added = states[t + 1] - p_minus
            worst["mass_err"] = max(worst["mass_err"],
                                    abs(float(states[t + 1].sum()) - 1.0))
            worst["min_bin"] = min(worst["min_bin"], float(states[t + 1].min()))
            worst["add_const_err"] = max(worst["add_const_err"],
                                         float(np.abs(added - r["removed_mass"] / K).max()))
            worst["redep_err"] = max(worst["redep_err"],
                                     abs(float(added.sum()) - r["removed_mass"]))
        v["worst"] = worst

    # ---------------- Stage 1 gate: A/B integrity vs committed M1C ------------------
    print("\n=== Stage 1: A/B baseline integrity vs committed M1C results ===")
    m1c = {(r["variant"], int(r["K"])): r
           for r in read_csv_rows(REPO_ROOT / "experiments/m1c_kernel_ablation/results/variant_summary.csv")}
    ann = {int(r["K"]): r
           for r in read_csv_rows(REPO_ROOT / "experiments/m1b_2_micro_scaling/results/annihilation_events.csv")}
    integ_rows: list[dict] = []
    fields = ["t_exact_stop", "t_U_best", "U_density_best", "U_density_at_stop",
              "improvement_retention", "moved_mass_median", "moved_mass_total",
              "D_max_p5", "D_max_p50", "D_max_p95", "U_density_p5",
              "U_density_p50", "U_density_p95", "share_delta_U_positive",
              "share_D_max_increases", "near_mass_05_at_stop"]
    for variant in ["A", "B"]:
        for K in K_list:
            for f in fields:
                exp = str(cache[(variant, K)][f])
                got = str(m1c[(variant, K)][f])
                integ_rows.append({"check": f"{variant} K={K}: {f}",
                                   "expected": got, "got": exp,
                                   "pass": exp == got})
    for K in K_list:
        stop_a = runs[("A", K)]["stop"]
        integ_rows.append({"check": f"A K={K}: t_exact_stop vs M1B.2",
                           "expected": ann[K]["t_stop"], "got": str(stop_a),
                           "pass": str(stop_a) == ann[K]["t_stop"]})
    n_pass = sum(1 for r in integ_rows if r["pass"])
    print(f"  baseline integrity: {n_pass}/{len(integ_rows)} checks passed")
    if n_pass < len(integ_rows):
        for r in integ_rows:
            if not r["pass"]:
                print(f"  [FAIL] {r['check']}: expected {r['expected']}, got {r['got']}")
        write_csv(out_dir / "baseline_integrity.csv",
                  ["check", "expected", "got", "pass"], integ_rows)
        print("BASELINE DRIFT — stopping before any A/U/B interpretation "
              "(config.baseline_policy). Resolve the A/B divergence first.")
        return 1

    # ---------------- per-round diagnostics (all three variants) --------------------
    per_rows: list[dict] = []
    for variant in ["A", "U", "B"]:
        for K in K_list:
            v = runs[(variant, K)]
            run_rows = v["run"]["rows"]
            n = v["n"]
            states = v["states"]
            cum = 0.0
            for t in range(n):
                r = run_rows[t]
                F_r = float((q_unif if variant == "U" else q_near)
                            .cdf(np.array([r["a_t"]]))[0])
                cum += r["removed_mass"]
                prefix_change = float(states[t + 1][: r["j_t"]].sum()
                                      - states[t][: r["j_t"]].sum())
                row = {
                    "variant": variant, "K": K, "t": t,
                    "active": r["active"], "j_t": r["j_t"], "a_t": r["a_t"],
                    "D_max": r["D_max"], "D_margin": r["D_margin"],
                    "sweep_mass": r["prefix_mass"], "moved_mass": r["removed_mass"],
                    "cum_moved_mass": cum,
                    "F_respray_a": F_r,
                    "prefix_mass_change": prefix_change,
                    "near_mass_05": r["near_mass_05"],
                    "L1_state_step": v["L1"][t],
                    "U_density": K * r["U_L2"],
                    "delta_U_density": K * (r["U_L2"] - run_rows[t + 1]["U_L2"]),
                }
                if variant in ["A", "U"]:
                    row["net_export_mass"] = r["removed_mass"] * (1.0 - F_r)
                else:
                    row["net_export_mass"] = ""
                per_rows.append({k: (f"{v2:.12e}" if isinstance(v2, float) else v2)
                                 for k, v2 in row.items()})
    write_csv(out_dir / "per_round_diagnostics.csv",
              list(per_rows[0]), per_rows)

    # ---------------- summary CSV ----------------
    summary_csv = out_dir / "variant_summary.csv"
    summary_fields: list[str] = []
    for row in summary_rows:
        for key in row:
            if key not in summary_fields:
                summary_fields.append(key)
    write_csv(summary_csv, summary_fields, summary_rows)
    write_csv(out_dir / "baseline_integrity.csv",
              ["check", "expected", "got", "pass"], integ_rows)
    print("\n=== variant summary ===")
    for r in summary_rows:
        print(f"  {r['variant']} K={r['K']:>3}: t_stop={r['t_exact_stop']:>5} "
              f"t_best={r['t_U_best']:>5}  U_stop={r['U_density_at_stop']}  "
              f"U_best={r['U_density_best']}  moved_tot={r['moved_mass_total']}  "
              f"Dmax_up={r['share_D_max_increases']}")
    print("\n=== U invariants ===")
    for K in K_list:
        w = runs[("U", K)]["worst"]
        print(f"  K={K}: |sum-1|<={w['mass_err']:.1e}  min_bin={w['min_bin']:.2e}  "
              f"|added-M/K|<={w['add_const_err']:.1e}  redep_err<={w['redep_err']:.1e}")

    # ---------------- epsilon overlay (all three variants) ---------------------------
    over_rows: list[dict] = []
    for variant in ["A", "U", "B"]:
        for K in K_list:
            v = runs[(variant, K)]
            n, Ud = v["n"], v["Ud"]
            run_rows = v["run"]["rows"]
            Dm_full = np.array([r["D_max"] for r in run_rows[:n + 1]])
            u0 = float(Ud[0])
            u_best = float(Ud[:n + 1].min())
            for eps in eps_list:
                hit = np.nonzero(Dm_full <= eps)[0]
                if not len(hit):
                    over_rows.append({"variant": variant, "K": K, "epsilon": eps,
                                      "t_stop": "none", "moved_mass_at_stop": "",
                                      "U_density_at_stop": "",
                                      "improvement_retention": ""})
                    continue
                t_eps = int(hit[0])
                moved_at = alpha * run_rows[t_eps]["prefix_mass"]
                over_rows.append({
                    "variant": variant, "K": K, "epsilon": eps,
                    "t_stop": t_eps,
                    "moved_mass_at_stop": f"{moved_at:.6f}",
                    "U_density_at_stop": f"{Ud[t_eps]:.6e}",
                    "improvement_retention":
                        f"{1.0 - (Ud[t_eps] - u_best) / (u0 - u_best):.6f}",
                })
    write_csv(out_dir / "epsilon_overlay.csv", list(over_rows[0]), over_rows)

    # ---------------- figures ----------------
    fig1 = out_dir / "fig_AUB_U_density_K100.png"
    plot_aub_trajectories(runs, K_primary=cfg["primary_K"], path=fig1)
    fig2 = out_dir / "fig_AUB_D_max_moved_K100.png"
    plot_aub_kernel_windows(runs, K_primary=cfg["primary_K"], path=fig2)

    # ---------------- metadata ----------------
    output_files = sorted(
        p.name for p in out_dir.iterdir() if p.is_file() and p.name != "metadata.json"
    )
    metadata = {
        "model": cfg["model"],
        "model_version": cfg["model_version"],
        "experiment": cfg["experiment"],
        "run_date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": cfg,
        "headline": {
            "baseline_integrity_checks_passed": f"{n_pass}/{len(integ_rows)}",
            "exact_stops": {variant: {str(K): runs[(variant, K)]["stop"]
                                      for K in K_list}
                            for variant in ["A", "U", "B"]},
            "U_stop_state_U_density": {
                str(K): float(runs[("U", K)]["Ud"][runs[("U", K)]["n"]])
                for K in K_list},
            "U_invariants_held": {str(K): runs[("U", K)]["worst"] for K in K_list},
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m1c1_uniform_kernel_ablation/run_experiment.py",
        "outputs": {name: sha256_of(out_dir / name) for name in output_files},
    }
    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"\noutputs written to: {out_dir} ({len(output_files) + 1} files)")
    print(f"runtime: {metadata['runtime_seconds']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
