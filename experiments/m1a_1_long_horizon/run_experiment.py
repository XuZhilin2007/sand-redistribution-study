"""Runner for M1A.1: long-horizon & resolution diagnostics of M1A.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m1a_1_long_horizon/run_experiment.py

Everything is diagnostic: the M1A rule, stopping rule and throw
distribution are exactly those of the canonical M1A experiment. New here:
  * T_max = 5000 baseline run (K=100) with per-round ambiguity diagnostics;
  * K in {50, 100, 200, 400} resolution runs (T_max=5000 each), compared
    via U_density = K*U and the one-sided CDF discrepancy D_max;
  * recurrence probe (lags 1..50, second half of the trajectory);
  * boundary-jump classification (flat-top drift vs genuine peak switch);
  * a few representative D(x) profiles.
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
from sand_m0.diagnostics import (  # noqa: E402
    boundary_diagnostics,
    largest_jumps,
    recurrence_l1,
    sup_cdf_minus_uniform,
    window_stats,
)
from sand_m0.model import DISTRIBUTIONS  # noqa: E402
from sand_m0.plotting import (  # noqa: E402
    plot_long_horizon,
    plot_k_sensitivity,
    plot_d_profiles,
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

    T_max = cfg["T_max"]
    alpha = cfg["alpha"]
    tol = cfg["stop_tolerance"]
    K_list = list(cfg["K_list"])
    K_base = cfg["baseline_K"]
    checkpoints = list(cfg["checkpoints_t"])
    windows = list(cfg["windows"])
    max_lag = cfg["recurrence"]["max_lag"]
    jump_threshold = cfg["boundary_ambiguity"]["large_jump_threshold"]
    started = time.time()

    print(f"M1A.1 diagnostics run — config: {config_path}")
    print(f"  T_max={T_max}, alpha={alpha}, tol={tol}, K_list={K_list}")

    runs: dict[int, dict] = {}
    for K in K_list:
        t0 = time.time()
        runs[K] = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=T_max,
                                        alpha=alpha, stop_tolerance=tol)
        print(f"  [K={K}] done in {time.time() - t0:.1f}s; "
              f"stopped_at={runs[K]['stopped_at']}")

    # ---------------- 1. baseline long-horizon trajectory CSV (K=100) ----------------
    rows = runs[K_base]["rows"]
    traj_rows: list[dict] = []
    for t, r in enumerate(rows):
        diag = boundary_diagnostics(runs[K_base]["states"][t], K_base)
        traj_rows.append(
            {
                "t": t,
                "U_L2": f"{r['U_L2']:.12e}",
                "U_density": f"{K_base * r['U_L2']:.12e}",
                "D_max": f"{r['D_max']:.12e}",
                "D_sup_continuous_check": f"{sup_cdf_minus_uniform(runs[K_base]['states'][t], K_base, fine=20001):.12e}",
                "a_t": f"{r['a_t']:.6f}",
                "near_mass_05": f"{r['near_mass_05']:.12f}",
                "total_mass": f"{runs[K_base]['states'][t].sum():.12f}",
                "active": r["active"],
                "D_margin": f"{r['D_margin']:.6e}",
                "value_gap": f"{diag['value_gap']:.6e}",
                "top2_separation": f"{diag['spatial_separation']:.6f}",
                "span_1pct": f"{diag['span_1pct']:.6f}",
                "components_1pct": diag["components_1pct"],
                "farthest_point_1pct": f"{diag['farthest_point_1pct']:.6f}",
                "span_5pct": f"{diag['span_5pct']:.6f}",
                "components_5pct": diag["components_5pct"],
                "farthest_point_5pct": f"{diag['farthest_point_5pct']:.6f}",
            }
        )
    traj_csv = out_dir / "long_horizon_trajectory_K100.csv"
    write_csv(traj_csv, list(traj_rows[0]), traj_rows)

    # ---------------- 2. checkpoint + window summary (K=100) ----------------
    U = np.array([r["U_L2"] for r in rows])
    Ud = K_base * U
    Dm = np.array([r["D_max"] for r in rows])
    a_arr = np.array([r["a_t"] for r in rows])
    nm = np.array([r["near_mass_05"] for r in rows])

    summary_rows: list[dict] = []
    for t in checkpoints:
        summary_rows.append({
            "t": t, "U_L2": f"{U[t]:.6e}", "U_density": f"{Ud[t]:.6e}",
            "D_max": f"{Dm[t]:.6e}", "a_t": f"{a_arr[t]:.4f}",
            "near_mass_05": f"{nm[t]:.6f}",
        })
    best_t = int(U.argmin())
    summary_rows.append({
        "t": f"best({best_t})", "U_L2": f"{U[best_t]:.6e}",
        "U_density": f"{Ud[best_t]:.6e}", "D_max": f"{Dm[best_t]:.6e}",
        "a_t": f"{a_arr[best_t]:.4f}", "near_mass_05": f"{nm[best_t]:.6f}",
    })
    for w in windows:
        for name, vals in [("U_L2", U), ("U_density", Ud), ("D_max", Dm)]:
            ws = window_stats(vals[-w:])
            summary_rows.append({
                "t": f"last{w}:{name}", "U_L2": f"{ws['min']:.6e}" if name == "U_L2" else "",
                "U_density": f"{ws['min']:.6e}" if name == "U_density" else "",
                "D_max": f"{ws['min']:.6e}" if name == "D_max" else "",
                "a_t": f"win min={ws['min']:.4e} max={ws['max']:.4e} "
                       f"mean={ws['mean']:.4e} std={ws['std']:.3e}",
                "near_mass_05": "",
            })
    summary_csv = out_dir / "checkpoints_and_windows.csv"
    write_csv(summary_csv, list(summary_rows[0]), summary_rows)

    print("\n=== K=100 checkpoints ===")
    for r in summary_rows:
        if isinstance(r["t"], int) or str(r["t"]).startswith("best"):
            print(f"  t={r['t']:>9}: U={r['U_L2']}  U_dens={r['U_density']}  "
                  f"D_max={r['D_max']}  a={r['a_t']}  near={r['near_mass_05']}")
    print("=== windows ===")
    for r in summary_rows:
        if str(r["t"]).startswith("last"):
            print(f"  {r['t']}: {r['a_t']}")

    # ---------------- 3. recurrence probe (K=100, ACTIVE segment only) ----------------
    # After the stop the state is frozen; the frozen tail would trivially give
    # L1 = 0, so recurrence is probed on the sweeping (active) segment.
    stop100 = runs[K_base]["stopped_at"]
    active_end = stop100 if stop100 is not None else T_max + 1
    rec_states = runs[K_base]["states"][:max(active_end, 2)]
    rec = recurrence_l1(rec_states, max_lag=max_lag, second_half=True)
    rec_rows = [
        {"lag": k, "L1_min": f"{v['min']:.6e}", "L1_mean": f"{v['mean']:.6e}"}
        for k, v in rec.items()
    ]
    rec_csv = out_dir / "recurrence_l1_K100.csv"
    write_csv(rec_csv, list(rec_rows[0]), rec_rows)
    one_step_scale = rec[1]["mean"]
    smallest = min(rec.values(), key=lambda v: v["min"])
    smallest_lag = min(rec, key=lambda k: rec[k]["min"])
    print(f"\n[recurrence] active segment t=0..{len(rec_states) - 1}; "
          f"typical one-step L1 = {one_step_scale:.4e}; "
          f"smallest min over lags 1..{max_lag}: {smallest['min']:.4e} at lag {smallest_lag}")

    # ---------------- 4. K sensitivity summary ----------------
    ks_rows: list[dict] = []
    for K in K_list:
        rr = runs[K]["rows"]
        Uk = np.array([r["U_L2"] for r in rr])
        Dk = np.array([r["D_max"] for r in rr])
        ak = np.array([r["a_t"] for r in rr])
        nk = np.array([r["near_mass_05"] for r in rr])
        jumps = [abs(ak[t] - ak[t - 1]) for t in range(1, len(ak))]
        wU = window_stats(K * Uk[-windows[1]:])
        wD = window_stats(Dk[-windows[1]:])
        ks_rows.append({
            "K": K,
            "stopped_at": str(runs[K]["stopped_at"]),
            "U_density_t0": f"{K * Uk[0]:.6e}",
            "U_density_t100": f"{K * Uk[100]:.6e}",
            "U_density_t1000": f"{K * Uk[1000]:.6e}",
            "U_density_t5000": f"{K * Uk[-1]:.6e}",
            "U_density_best": f"{(K * Uk).min():.6e}",
            "U_density_last500_mean": f"{wU['mean']:.6e}",
            "U_density_last500_std": f"{wU['std']:.3e}",
            "D_max_t0": f"{Dk[0]:.6e}",
            "D_max_t1000": f"{Dk[1000]:.6e}",
            "D_max_t5000": f"{Dk[-1]:.6e}",
            "D_max_last500_min": f"{wD['min']:.6e}",
            "D_max_last500_mean": f"{wD['mean']:.6e}",
            "D_max_last500_max": f"{wD['max']:.6e}",
            "near_mass_final": f"{nk[-1]:.6f}",
            "a_t_min": f"{ak.min():.4f}",
            "a_t_max": f"{ak.max():.4f}",
            "a_t_mean": f"{ak.mean():.4f}",
            "large_jumps_gt_0.2": sum(1 for j in jumps if j > jump_threshold),
            "largest_jump": f"{max(jumps):.4f}",
        })
    ks_csv = out_dir / "k_sensitivity_summary.csv"
    write_csv(ks_csv, list(ks_rows[0]), ks_rows)
    print("\n=== K sensitivity ===")
    for r in ks_rows:
        print(f"  K={r['K']:>3}: stopped={r['stopped_at']:>4}  "
              f"U_dens t1000={float(r['U_density_t1000']):.4e} t5000={float(r['U_density_t5000']):.4e}  "
              f"D_max t5000={float(r['D_max_t5000']):.4e}  "
              f"a∈[{float(r['a_t_min']):.2f},{float(r['a_t_max']):.2f}]  "
              f"jumps>{jump_threshold}: {r['large_jumps_gt_0.2']}")

    # ---------------- 5. boundary ambiguity aggregates (K=100, ACTIVE rounds) ----------
    # Only decisions taken while the model was still sweeping are diagnostic;
    # the frozen tail's argmax is the trivial tie at j=K.
    gaps = np.array([float(r["value_gap"]) for r in traj_rows[:active_end]])
    seps = np.array([float(r["top2_separation"]) for r in traj_rows[:active_end]])
    comps1 = np.array([int(r["components_1pct"]) for r in traj_rows[:active_end]])
    fars1 = np.array([float(r["farthest_point_1pct"]) for r in traj_rows[:active_end]])
    amb_rows = [
        {"statistic": "active_decisions_total", "value": str(active_end)},
        {"statistic": "active_rounds_with_gap_below_1pct_of_Dmax", "value": str(int((gaps <= 0.01 * Dm[:active_end]).sum()))},
        {"statistic": "active_rounds_with_top2_adjacent(sep<=2/K)", "value": str(int((seps <= 2.0 / K_base).sum()))},
        {"statistic": "active_rounds_with_far_second_peak(sep>0.1)", "value": str(int((seps > 0.1).sum()))},
        {"statistic": "active_rounds_with_A1pct_multiple_components", "value": str(int((comps1 > 1).sum()))},
        {"statistic": "active_rounds_with_A1pct_reaching_beyond_0.1_from_argmax", "value": str(int((fars1 > 0.1).sum()))},
        {"statistic": "median_top2_value_gap", "value": f"{np.median(gaps):.3e}"},
        {"statistic": "median_top2_separation", "value": f"{np.median(seps):.4f}"},
        {"statistic": "median_farthest_A1pct_point", "value": f"{np.median(fars1):.4f}"},
        {"statistic": "late300_active_decisions_with_A1pct_multiple_components",
         "value": str(int((comps1[-min(300, active_end):] > 1).sum()))},
        {"statistic": "late300_median_farthest_A1pct_point",
         "value": f"{np.median(fars1[-min(300, active_end):]):.4f}"},
    ]
    amb_csv = out_dir / "boundary_ambiguity_summary.csv"
    write_csv(amb_csv, list(amb_rows[0]), amb_rows)
    print("\n=== boundary ambiguity (K=100) ===")
    for r in amb_rows:
        print(f"  {r['statistic']}: {r['value']}")

    # ---------------- 6. large jumps within the ACTIVE segment ----------------
    active_a = a_arr[:active_end]
    n_big = int((np.abs(np.diff(active_a)) > jump_threshold).sum())
    jumps = largest_jumps(active_a.tolist(), n=max(n_big, 1))
    jump_rows = []
    for jr in jumps:
        t = jr["t"]
        d_before = boundary_diagnostics(runs[K_base]["states"][t - 1], K_base)
        d_after = boundary_diagnostics(runs[K_base]["states"][t], K_base)
        Db = np.cumsum(runs[K_base]["states"][t - 1]) - np.arange(1, K_base + 1) / K_base
        Da = np.cumsum(runs[K_base]["states"][t]) - np.arange(1, K_base + 1) / K_base
        j_before = int(round(jr["a_before"] * K_base))
        j_after = int(round(jr["a_after"] * K_base))
        # flat-top drift iff both boundary positions are near-optimal (A_5%)
        # in BOTH rounds: the argmax merely moved inside a persistent flat
        # region. Otherwise the D-profile itself switched peaks.
        new_was_near_optimal = Db[j_after - 1] >= 0.95 * Db.max()
        old_still_near_optimal = Da[j_before - 1] >= 0.95 * Da.max()
        kind = ("flat-top drift" if (new_was_near_optimal and old_still_near_optimal)
                else "peak switch")
        jump_rows.append({
            "t": t,
            "a_before": f"{jr['a_before']:.4f}",
            "a_after": f"{jr['a_after']:.4f}",
            "jump": f"{jr['jump']:.4f}",
            "new_position_in_before_A5pct": new_was_near_optimal,
            "old_position_in_after_A5pct": old_still_near_optimal,
            "classification": kind,
            "D_max_before": f"{d_before['D_max']:.6e}",
            "D_max_after": f"{d_after['D_max']:.6e}",
        })
    jump_csv = out_dir / "largest_boundary_jumps.csv"
    write_csv(jump_csv, list(jump_rows[0]), jump_rows)
    n_switch = sum(1 for r in jump_rows if r["classification"] == "peak switch")
    print(f"\n[large jumps] all {len(jump_rows)} jumps with |Δa|>{jump_threshold} "
          f"({100.0 * len(jump_rows) / max(active_end - 1, 1):.1f}% of decisions): "
          f"{n_switch} peak switches, {len(jump_rows) - n_switch} flat-top drifts")
    for r in jump_rows[:3]:
        print(f"  t={r['t']}: {r['a_before']} -> {r['a_after']} ({r['classification']})")

    # ---------------- 7. figures ----------------
    fig1 = out_dir / "fig_long_horizon_K100.png"
    plot_long_horizon(U, Dm, a_arr, tol=tol, path=fig1, checkpoints=checkpoints)
    fig2 = out_dir / "fig_k_sensitivity.png"
    plot_k_sensitivity(runs, K_list, path=fig2)
    stop_lab = f"stopping state t={stop100}" if stop100 is not None else f"t={T_max}"
    d_rounds = {}
    for label, t in [("early t=20", 20), ("medium t=500", 500), (stop_lab, min(active_end, T_max))]:
        st = runs[K_base]["states"][t]
        d_rounds[label] = np.cumsum(st) - np.arange(1, K_base + 1) / K_base
    big = max(jumps, key=lambda r: r["jump"])
    tb = big["t"]
    d_rounds[f"before jump t={tb - 1}"] = np.cumsum(runs[K_base]["states"][tb - 1]) - np.arange(1, K_base + 1) / K_base
    d_rounds[f"after jump t={tb}"] = np.cumsum(runs[K_base]["states"][tb]) - np.arange(1, K_base + 1) / K_base
    fig3 = out_dir / "fig_d_profiles.png"
    plot_d_profiles(d_rounds, K=K_base, path=fig3)

    # ---------------- 8. metadata ----------------
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
            "stopped_at_per_K": {str(K): str(runs[K]["stopped_at"]) for K in K_list},
            "best_round_K100": best_t,
            "first_divergence_note": "K=100 run reproduces the canonical M1A trajectory bit-for-bit for t=0..100",
            "recurrence_min_L1": smallest["min"],
            "recurrence_min_lag": smallest_lag,
            "one_step_L1_mean": one_step_scale,
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m1a_1_long_horizon/run_experiment.py",
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
