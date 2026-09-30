"""Runner for M1C: correction-kernel ablation (Variant A vs Variant B).

One command regenerates every CSV, figure and the metadata file:

    python experiments/m1c_kernel_ablation/run_experiment.py

Research question (frozen): does churn persist when the M1 selection rule is
held fixed but the fixed biased redistribution law is replaced by a
target-aware corrective kernel?

  * Variant A — canonical M1A/M1B baseline, untouched: removed mass is
    resprayed with the fixed throw distribution q_near.
  * Variant B — selection identical (same cumulative-excess argmax boundary,
    same alpha removal, same stopping rules); the removed mass is
    redistributed by target-aware deficit filling against u_i = 1/K.

Stage 1 of the runner REPRODUCES the committed M1B.2/M1A.1 baseline with
Variant A and hard-stops if it diverges; only then are A-vs-B comparisons
produced. No post-hoc tuning: each variant runs once per K under the frozen
configuration and the trajectories are accepted as-is.
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
from sand_m0.plotting import plot_a_vs_b_trajectories, plot_a_vs_b_summary  # noqa: E402


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
    q = DISTRIBUTIONS["near"]
    started = time.time()

    print(f"M1C correction-kernel ablation — config: {config_path}")
    print(f"  K={K_list}, alpha={alpha}, exact tol={tol}, T_max={T_max}")

    # ---------------- 0. Stage 1: baseline reproduction with Variant A ----------------
    runs: dict[tuple[str, int], dict] = {}
    print("\n=== Stage 1: Variant A baseline reproduction ===")
    base_ok = True
    repro_rows: list[dict] = []

    def record(check: str, expected: float | int | str, got: float | int | str,
               ok: bool, note: str = "") -> None:
        nonlocal base_ok
        base_ok = base_ok and ok
        repro_rows.append({
            "check": check, "expected": expected, "got": got,
            "pass": ok, **({"note": note} if note else {}),
        })
        flag = "OK " if ok else "FAIL"
        print(f"  [{flag}] {check}: expected {expected}, got {got}")

    m1b2 = REPO_ROOT / "experiments/m1b_2_micro_scaling/results"
    ann = {int(r["K"]): r for r in read_csv_rows(m1b2 / "annihilation_events.csv")}
    dipstat = {int(r["K"]): r for r in read_csv_rows(m1b2 / "band_and_dip_statistics.csv")}
    eps_overlay_ref = {(int(r["K"]), float(r["epsilon"])): r
                       for r in read_csv_rows(m1b2 / "epsilon_overlay.csv")}
    scale_ref = {r["quantity"]: r for r in read_csv_rows(m1b2 / "cross_K_scaling.csv")}

    for K in K_list:
        t0 = time.time()
        run = run_m1a_deterministic(q, K=K, T_max=T_max, alpha=alpha,
                                    stop_tolerance=tol, redistribution="throw")
        runs[("A", K)] = {"run": run, "stop": run["stopped_at"]}
        stop = run["stopped_at"]
        rows = run["rows"]
        states = np.asarray(run["states"], dtype=float)
        Ud = np.array([K * r["U_L2"] for r in rows])
        Dm = np.array([r["D_max"] for r in rows[:stop]])
        moved = np.array([r["removed_mass"] for r in rows[:stop]])
        F_q_a = np.array([float(q.cdf(np.array([r["a_t"]]))[0]) for r in rows[:stop]])
        L1 = np.abs(states[1:stop + 1] - states[:stop]).sum(axis=1)
        net_export = moved * (1.0 - F_q_a)
        runs[("A", K)].update({"Ud": Ud, "Dm": Dm, "moved": moved, "L1": L1,
                               "net_export": net_export})

        # exact stopping + annihilation structure vs committed values
        record(f"A K={K}: t_exact_stop", ann[K]["t_stop"], stop, str(stop) == ann[K]["t_stop"])
        pre = rows[stop - 1]
        record(f"A K={K}: pre-stop a_t (4dp)", ann[K]["pre_stop_a_t"], f"{pre['a_t']:.4f}",
               f"{pre['a_t']:.4f}" == ann[K]["pre_stop_a_t"])
        record(f"A K={K}: pre-stop moved (6dp)", ann[K]["moved_mass"], f"{moved[stop-1]:.6f}",
               f"{moved[stop-1]:.6f}" == ann[K]["moved_mass"])

        # late-stage (D_max < 0.005) medians vs committed cross_K_scaling.csv (6dp)
        late = Dm < 0.005
        for name, arr in [("moved_mass", moved), ("net_export_mass", net_export),
                          ("L1_state_step", L1)]:
            med = float(np.median(arr[late])) if late.any() else float("nan")
            record(f"A K={K}: late median {name} (6dp)", scale_ref[name][f"late_median_K{K}"],
                   f"{med:.6f}", f"{med:.6f}" == scale_ref[name][f"late_median_K{K}"])

        # fluctuation bands + dips vs committed band_and_dip_statistics.csv (5dp)
        ref = dipstat[K]
        for name, arr in [("D_max", Dm), ("U_density", Ud[:stop])]:
            for p in [5, 50, 95]:
                key = f"{name}_p{p}"
                val = f"{float(np.percentile(arr, p)):.5f}"
                record(f"A K={K}: {key} (5dp)", ref[key], val, val == ref[key])
        for eps in eps_list:
            dips = int((Dm <= eps).sum())
            first = int(np.nonzero(Dm <= eps)[0][0]) if dips else -1
            record(f"A K={K}: dips<= {eps}", ref[f"dips_le_{eps}"], dips,
                   dips == int(ref[f"dips_le_{eps}"]))
            record(f"A K={K}: first_dip {eps}", ref[f"first_dip_{eps}"], first,
                   first == int(ref[f"first_dip_{eps}"]))

        # epsilon first-passage vs committed epsilon_overlay.csv
        u_best = float(Ud[:stop].min())
        for eps in eps_list:
            hit = np.nonzero(Dm <= eps)[0]
            if not len(hit):
                continue
            t_eps = int(hit[0])
            ref_r = eps_overlay_ref[(K, eps)]
            record(f"A K={K} eps={eps}: t_stop", ref_r["t_stop"], t_eps,
                   str(t_eps) == ref_r["t_stop"])
            ret = 1.0 - (Ud[t_eps] - u_best) / (Ud[0] - u_best)
            record(f"A K={K} eps={eps}: retention (6dp)",
                   ref_r["improvement_retention_fine"], f"{ret:.6f}",
                   f"{ret:.6f}" == ref_r["improvement_retention_fine"])

        print(f"  [K={K}] A reproduced ({time.time() - t0:.1f}s)")

    repro_csv = out_dir / "baseline_reproduction.csv"
    write_csv(repro_csv, ["check", "expected", "got", "pass", "note"],
              [{k: str(v) for k, v in r.items()} for r in repro_rows])
    n_checks = len(repro_rows)
    n_pass = sum(1 for r in repro_rows if r["pass"])
    print(f"  baseline reproduction: {n_pass}/{n_checks} checks passed")
    if not base_ok:
        print("BASELINE DIVERGENCE — stopping before any A-vs-B interpretation "
              "(config.baseline_policy). Resolve the Variant A divergence first.")
        return 1

    # ---------------- 1. Variant B runs + per-round diagnostics ----------------
    print("\n=== Stage 2: Variant B (deficit-fill kernel) ===")
    for K in K_list:
        t0 = time.time()
        run = run_m1a_deterministic(q, K=K, T_max=T_max, alpha=alpha,
                                    stop_tolerance=tol, redistribution="deficit_fill")
        stop = run["stopped_at"]
        rows = run["rows"]
        states = np.asarray(run["states"], dtype=float)
        Ud = np.array([K * r["U_L2"] for r in rows])
        Dm = np.array([r["D_max"] for r in rows[:stop]])
        moved = np.array([r["removed_mass"] for r in rows[:stop]])
        L1 = np.abs(states[1:stop + 1] - states[:stop]).sum(axis=1)
        runs[("B", K)] = {"run": run, "stop": stop, "Ud": Ud, "Dm": Dm,
                          "moved": moved, "L1": L1, "states": states}

        # round-by-round kernel diagnostics + invariant verification
        worst = {"mass_err": 0.0, "min_bin": np.inf, "overfill": -np.inf,
                 "identity_err": 0.0, "fill_ratio": -np.inf}
        for t in range(stop):
            r = rows[t]
            j_t = r["j_t"]
            p_minus = states[t].copy()
            p_minus[:j_t] *= (1.0 - alpha)
            p_next = states[t + 1]
            deficit = np.maximum(1.0 / K - p_minus, 0.0)
            excess = np.maximum(p_minus - 1.0 / K, 0.0)
            r["total_deficit"] = float(deficit.sum())
            r["post_removal_excess"] = float(excess.sum())
            r["fill_ratio"] = r["removed_mass"] / r["total_deficit"]
            worst["mass_err"] = max(worst["mass_err"],
                                    abs(float(states[t + 1].sum()) - 1.0))
            worst["min_bin"] = min(worst["min_bin"], float(p_next.min()))
            worst["overfill"] = max(worst["overfill"],
                                    float((p_next - 1.0 / K).max()))
            worst["identity_err"] = max(
                worst["identity_err"],
                abs(r["total_deficit"] - r["post_removal_excess"] - r["removed_mass"]))
            worst["fill_ratio"] = max(worst["fill_ratio"], r["fill_ratio"])
        runs[("B", K)]["worst"] = worst
        print(f"  [K={K}] B stop t={stop}; invariants: |sum-1|<={worst['mass_err']:.1e}, "
              f"min_bin={worst['min_bin']:.3e}, overfill<={worst['overfill']:.1e}, "
              f"identity_err<={worst['identity_err']:.1e}, "
              f"fill_ratio<={worst['fill_ratio']:.6f} ({time.time() - t0:.1f}s)")

    # ---------------- 2. per-round diagnostics CSVs (both variants) ----------------
    per_rows: list[dict] = []
    for variant in ["A", "B"]:
        for K in K_list:
            v = runs[(variant, K)]
            run_rows = v["run"]["rows"]
            stop = v["stop"]
            states = np.asarray(v["run"]["states"], dtype=float)
            for t in range(stop):
                r = run_rows[t]
                F_q_a = float(q.cdf(np.array([r["a_t"]]))[0])
                prefix_change = float(states[t + 1][:r["j_t"]].sum()
                                      - states[t][:r["j_t"]].sum())
                row = {
                    "variant": variant, "K": K, "t": t,
                    "active": r["active"], "j_t": r["j_t"], "a_t": r["a_t"],
                    "D_max": r["D_max"], "D_margin": r["D_margin"],
                    "sweep_mass": r["prefix_mass"], "moved_mass": r["removed_mass"],
                    "F_q_a": F_q_a,
                    "prefix_mass_change": prefix_change,
                    "near_mass_05": r["near_mass_05"],
                    "L1_state_step": v["L1"][t],
                    "U_density": K * r["U_L2"],
                    "delta_U_density": K * (r["U_L2"] - run_rows[t + 1]["U_L2"]),
                }
                if variant == "A":
                    row["net_export_mass"] = r["removed_mass"] * (1.0 - F_q_a)
                    row["total_deficit"] = ""
                    row["post_removal_excess"] = ""
                    row["fill_ratio"] = ""
                else:
                    row["net_export_mass"] = ""
                    row["total_deficit"] = r["total_deficit"]
                    row["post_removal_excess"] = r["post_removal_excess"]
                    row["fill_ratio"] = r["fill_ratio"]
                per_rows.append({k: (f"{v2:.12e}" if isinstance(v2, float) else v2)
                                 for k, v2 in row.items()})
    per_csv = out_dir / "per_round_diagnostics.csv"
    write_csv(per_csv, list(per_rows[0]), per_rows)

    # ---------------- 3. variant summary (churn / convergence quantities) ----------------
    def band(v: np.ndarray) -> dict:
        return {"p5": float(np.percentile(v, 5)), "p50": float(np.percentile(v, 50)),
                "p95": float(np.percentile(v, 95))}

    summary_rows: list[dict] = []
    for variant in ["A", "B"]:
        for K in K_list:
            v = runs[(variant, K)]
            stop, Ud, Dm, moved = v["stop"], v["Ud"], v["Dm"], v["moved"]
            # the stopping state is state `stop` (rows are state-indexed: row t
            # holds state t's metrics and the decision for t -> t+1)
            traj_Ud = Ud[:stop + 1]
            u0, u_best = float(traj_Ud[0]), float(traj_Ud.min())
            t_best = int(traj_Ud.argmin())
            u_stop = float(traj_Ud[-1])
            Db = band(Dm)
            Ub = band(traj_Ud)
            n = stop
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
                "share_D_max_increases":
                    f"{float((np.diff(Dm) > 0).mean()):.4f}",
                "near_mass_05_at_stop":
                    f"{v['run']['rows'][stop]['near_mass_05']:.6f}",
            }
            if variant == "B":
                w = v["worst"]
                row.update({
                    "inv_mass_err_max": f"{w['mass_err']:.3e}",
                    "inv_min_bin": f"{w['min_bin']:.3e}",
                    "inv_overfill_max": f"{w['overfill']:.3e}",
                    "inv_identity_err_max": f"{w['identity_err']:.3e}",
                    "inv_fill_ratio_max": f"{w['fill_ratio']:.6f}",
                })
            summary_rows.append(row)
    summary_csv = out_dir / "variant_summary.csv"
    summary_fields: list[str] = []
    for row in summary_rows:
        for key in row:
            if key not in summary_fields:
                summary_fields.append(key)
    write_csv(summary_csv, summary_fields, summary_rows)
    print("\n=== variant summary ===")
    for r in summary_rows:
        print(f"  {r['variant']} K={r['K']:>3}: t_stop={r['t_exact_stop']:>5} "
              f"t_best={r['t_U_best']:>5}  U_stop={r['U_density_at_stop']}  "
              f"retention={r['improvement_retention']}  "
              f"moved_median={r['moved_mass_median']}  "
              f"share_dU>0={r['share_delta_U_positive']}")

    # ---------------- 4. epsilon overlay (both variants, first passage) ----------------
    over_rows: list[dict] = []
    for variant in ["A", "B"]:
        for K in K_list:
            v = runs[(variant, K)]
            stop, Ud = v["stop"], v["Ud"]
            run_rows = v["run"]["rows"]
            # first passage over the FULL stopped trajectory (states 0..stop,
            # including the stopping state itself)
            Dm_full = np.array([r["D_max"] for r in run_rows[:stop + 1]])
            u0 = float(Ud[0])
            u_best = float(Ud[:stop + 1].min())
            for eps in eps_list:
                hit = np.nonzero(Dm_full <= eps)[0]
                if not len(hit):
                    over_rows.append({"variant": variant, "K": K, "epsilon": eps,
                                      "t_stop": "none", "moved_mass_at_stop": "",
                                      "U_density_at_stop": "",
                                      "improvement_retention": ""})
                    continue
                t_eps = int(hit[0])
                # gross action scale of the decision at the stopping state
                # (= removed_mass for rounds that executed, M1B.2 convention)
                moved_at = alpha * run_rows[t_eps]["prefix_mass"]
                over_rows.append({
                    "variant": variant, "K": K, "epsilon": eps,
                    "t_stop": t_eps,
                    "moved_mass_at_stop": f"{moved_at:.6f}",
                    "U_density_at_stop": f"{Ud[t_eps]:.6e}",
                    "improvement_retention":
                        f"{1.0 - (Ud[t_eps] - u_best) / (u0 - u_best):.6f}",
                })
    over_csv = out_dir / "epsilon_overlay.csv"
    write_csv(over_csv, list(over_rows[0]), over_rows)

    # ---------------- 5. figures ----------------
    fig1 = out_dir / "fig_A_vs_B_K100.png"
    plot_a_vs_b_trajectories(runs, K_primary=cfg["primary_K"], path=fig1)
    fig2 = out_dir / "fig_A_vs_B_summary.png"
    plot_a_vs_b_summary(runs, K_list, path=fig2)

    # ---------------- 6. metadata ----------------
    output_files = sorted(
        p.name for p in out_dir.iterdir() if p.is_file() and p.name != "metadata.json"
    )
    a_stops = {K: runs[("A", K)]["stop"] for K in K_list}
    b_stops = {K: runs[("B", K)]["stop"] for K in K_list}
    metadata = {
        "model": cfg["model"],
        "model_version": cfg["model_version"],
        "experiment": cfg["experiment"],
        "run_date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": cfg,
        "headline": {
            "baseline_reproduction_checks_passed": f"{n_pass}/{n_checks}",
            "variant_A_exact_stops": a_stops,
            "variant_B_exact_stops": b_stops,
            "variant_B_post_stop_U_density": {
                str(K): float(runs[("B", K)]["Ud"][runs[("B", K)]["stop"]])
                for K in K_list},
            "variant_B_invariants_held": {
                str(K): runs[("B", K)]["worst"] for K in K_list},
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m1c_kernel_ablation/run_experiment.py",
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
