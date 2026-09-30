"""Runner for M1A: adaptive sweep boundary vs fixed-boundary M0.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m1a_adaptive_boundary/run_experiment.py

Contents (all parameters frozen in config.json):
  1. baseline M1A run (q_near, alpha=0.25, T_max=100) with boundary and
     excess trajectories;
  2. same-code M0 comparison runs (fixed a=0.5, T=100) for q_near and the
     two controls;
  3. controls: uniform and far-biased throw distributions under M1A
     (both are expected to stop immediately — no cumulative near excess);
  4. alpha structural check {0.10, 0.25, 0.50, 1.00}: boundary sequences,
     stopping, first divergence;
  5. lightweight MC spot check (5 seeds) including how often the noisy
     histogram picks a different boundary than the deterministic state.
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

from sand_m0.adaptive import (  # noqa: E402
    boundary_decision,
    cumulative_excess,
    run_m1a_deterministic,
    run_m1a_mc,
)
from sand_m0.model import DISTRIBUTIONS, noise_floor_l2  # noqa: E402
from sand_m0.plotting import (  # noqa: E402
    plot_initial_excess,
    plot_m0_vs_m1a,
    plot_boundary_and_excess,
)
from sand_m0.simulate import run_deterministic, run_deterministic_states  # noqa: E402


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


def trajectory_summary(rows: list[dict]) -> dict:
    best_t = min(range(len(rows)), key=lambda t: rows[t]["U_L2"])
    active_ts = [r["t"] for r in rows if r["active"]]
    return {
        "best_round": best_t,
        "best_U": rows[best_t]["U_L2"],
        "U_final": rows[-1]["U_L2"],
        "U_initial": rows[0]["U_L2"],
        "sweeps_executed": len(active_ts),
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

    K, T_max, N = cfg["K"], cfg["T_max"], cfg["N"]
    alpha0 = cfg["baseline_alpha"]
    tol = cfg["adaptive_rule"]["stop_tolerance"]
    alphas = list(cfg["alpha_structural_check"]["alphas"])
    seeds = list(cfg["mc_spotcheck"]["seeds"])
    started = time.time()

    print(f"M1A run — config: {config_path}")
    print(f"  adaptive rule: {cfg['adaptive_rule']['definition']}")
    print(f"  alpha={alpha0} (baseline), K={K}, T_max={T_max}, tol={tol}")

    # ---------------- 1. baseline M1A (q_near) ----------------
    m1a = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=T_max,
                                alpha=alpha0, stop_tolerance=tol)
    m1a_csv = out_dir / "m1a_baseline_trajectory.csv"
    write_csv(m1a_csv, list(m1a["rows"][0]), m1a["rows"])
    s = trajectory_summary(m1a["rows"])
    print(f"\n[M1A q_near] stopped_at={m1a['stopped_at']}  best t={s['best_round']} "
          f"U={s['best_U']:.6e}  U_final={s['U_final']:.6e}  "
          f"sweeps={s['sweeps_executed']}")
    a_arr = [r["a_t"] for r in m1a["rows"]]

    # ---------------- 2. M0 comparison (same code, T=T_max) ----------------
    m0 = run_deterministic(DISTRIBUTIONS["near"], K=K, T=T_max, a=0.5, alpha=alpha0)
    m0_u = [r["U_L2"] for r in m0]
    m0_near = [r["near_mass"] for r in m0]
    m0_best = min(range(T_max + 1), key=lambda t: m0_u[t])
    print(f"[M0  q_near] best t={m0_best} U={m0_u[m0_best]:.6e}  "
          f"U_final={m0_u[-1]:.6e}  near_mass_final={m0_near[-1]:.6f}")

    first_div = next((t for t in range(T_max + 1)
                      if abs(m1a["rows"][t]["U_L2"] - m0_u[t]) > 1e-15), None)
    print(f"first round where M1A U leaves the M0 trajectory: t={first_div}")

    comp_rows = [
        {
            "model": "M0 (fixed a=0.5)",
            "distribution": "q_near",
            "alpha": alpha0,
            "best_round": m0_best,
            "best_U": f"{m0_u[m0_best]:.12e}",
            "U_final": f"{m0_u[-1]:.12e}",
            "U_initial": f"{m0_u[0]:.12e}",
            "final_worse_than_best": m0_u[-1] > m0_u[m0_best],
            "near_mass_final": f"{m0_near[-1]:.10f}",
            "near_mass_min": f"{min(m0_near):.10f}",
            "sweeps_executed": T_max,
            "stopped_at": "never (fixed rule always sweeps)",
        },
        {
            "model": "M1A (adaptive)",
            "distribution": "q_near",
            "alpha": alpha0,
            "best_round": s["best_round"],
            "best_U": f"{s['best_U']:.12e}",
            "U_final": f"{s['U_final']:.12e}",
            "U_initial": f"{s['U_initial']:.12e}",
            "final_worse_than_best": s["U_final"] > s["best_U"],
            "near_mass_final": f"{m1a['rows'][-1]['near_mass_05']:.10f}",
            "near_mass_min": f"{min(r['near_mass_05'] for r in m1a['rows']):.10f}",
            "sweeps_executed": s["sweeps_executed"],
            "stopped_at": str(m1a["stopped_at"]),
        },
    ]

    # ---------------- 3. controls under M1A ----------------
    for name in cfg["controls"]:
        res = run_m1a_deterministic(DISTRIBUTIONS[name], K=K, T_max=10,
                                    alpha=alpha0, stop_tolerance=tol)
        r0 = res["rows"][0]
        comp_rows.append(
            {
                "model": "M1A (adaptive)",
                "distribution": f"q_{name}",
                "alpha": alpha0,
                "best_round": 0,
                "best_U": f"{r0['U_L2']:.12e}",
                "U_final": f"{res['rows'][-1]['U_L2']:.12e}",
                "U_initial": f"{r0['U_L2']:.12e}",
                "final_worse_than_best": False,
                "near_mass_final": f"{r0['near_mass_05']:.10f}",
                "near_mass_min": f"{r0['near_mass_05']:.10f}",
                "sweeps_executed": 0,
                "stopped_at": str(res["stopped_at"]),
            }
        )
        print(f"[M1A q_{name}] stopped_at={res['stopped_at']}  "
              f"D_max(0)={r0['D_max']:.2e}  sweeps=0")
    comp_csv = out_dir / "m0_vs_m1a_comparison.csv"
    write_csv(comp_csv, list(comp_rows[0]), comp_rows)

    # ---------------- 4. alpha structural check ----------------
    alpha_rows: list[dict] = []
    a_seqs: dict[float, list[float]] = {}
    u_seqs: dict[float, list[float]] = {}
    stops: dict[float, object] = {}
    sums: dict[float, dict] = {}
    for alpha in alphas:
        res = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=T_max,
                                    alpha=alpha, stop_tolerance=tol)
        rows = res["rows"]
        a_seqs[alpha] = [r["a_t"] for r in rows]
        u_seqs[alpha] = [r["U_L2"] for r in rows]
        stops[alpha] = res["stopped_at"]
        sums[alpha] = trajectory_summary(rows)
    for alpha in alphas:
        rows_a = a_seqs[alpha]
        sbest = sums[alpha]
        first_div = next(
            (t for t in range(1, T_max + 1)
             if abs(rows_a[t] - a_seqs[alpha0][t]) > 1e-12),
            None,
        )
        alpha_rows.append(
            {
                "alpha": alpha,
                "a_at_t0": rows_a[0],
                "a_at_t1": rows_a[1],
                "a_at_t2": rows_a[2],
                "first_divergence_round_vs_baseline": "" if first_div is None else first_div,
                "best_round": sbest["best_round"],
                "best_U": f"{sbest['best_U']:.12e}",
                "U_final": f"{sbest['U_final']:.12e}",
                "sweeps_executed": sbest["sweeps_executed"],
                "stopped_at": str(stops[alpha]),
                "a_t_min": min(rows_a),
                "a_t_max": max(rows_a),
            }
        )
        print(f"[alpha={alpha}] a(1)={rows_a[1]}  first_div={first_div}  "
              f"best t={sbest['best_round']} U={sbest['best_U']:.6e}  "
              f"stopped={stops[alpha]}  a range=[{min(rows_a):.2f},{max(rows_a):.2f}]")
    alpha_csv = out_dir / "alpha_structural_check.csv"
    write_csv(alpha_csv, list(alpha_rows[0]), alpha_rows)

    # state-level divergence: first round where alpha runs differ in U
    state_div = {
        alpha: next((t for t in range(T_max + 1)
                     if abs(u_seqs[alpha][t] - u_seqs[alpha0][t]) > 1e-15), None)
        for alpha in alphas if alpha != alpha0
    }
    print(f"first U divergence vs baseline: {state_div}")

    # ---------------- 5. MC spot check ----------------
    floor = noise_floor_l2(N, K)
    mc_rows: list[dict] = []
    boundary_flips = []
    det_a = [r["a_t"] for r in m1a["rows"]]
    for seed in seeds:
        res = run_m1a_mc(DISTRIBUTIONS["near"], K=K, T_max=T_max, alpha=alpha0,
                         stop_tolerance=tol, N=N, seed=seed)
        flips = sum(1 for r in res["rows"] if r["a_t"] != det_a[r["t"]])
        boundary_flips.append(flips)
        mc_rows.append(
            {
                "seed": seed,
                "U_final": f"{res['rows'][-1]['U_L2']:.12e}",
                "U_min": f"{min(r['U_L2'] for r in res['rows']):.12e}",
                "boundary_flips_vs_deterministic": flips,
            }
        )
    mc_csv = out_dir / "mc_spotcheck.csv"
    write_csv(mc_csv, list(mc_rows[0]), mc_rows)
    print(f"[MC] boundary flips vs deterministic per seed: {boundary_flips} "
          f"(out of {T_max} decisions); det final U={m1a['rows'][-1]['U_L2']:.6e}")

    # ---------------- 6. figures ----------------
    fig1 = out_dir / "fig_m0_vs_m1a_uniformity.png"
    plot_m0_vs_m1a(m0_u, [r["U_L2"] for r in m1a["rows"]], m0_best,
                   trajectory_summary(m1a["rows"])["best_round"], path=fig1)
    fig2 = out_dir / "fig_boundary_and_excess.png"
    plot_boundary_and_excess(a_seqs, u_seqs, {r["t"]: r["D_max"] for r in m1a["rows"]},
                             alpha0=alpha0, tol=tol, path=fig2)
    fig3 = out_dir / "fig_initial_cumulative_excess.png"
    plot_initial_excess(
        {name: cumulative_excess(DISTRIBUTIONS[name].bin_probs(K), K)
         for name in ["near", "uniform", "far"]},
        K=K, path=fig3,
    )

    # ---------------- 7. metadata ----------------
    output_files = sorted(
        p.name for p in out_dir.iterdir() if p.is_file() and p.name != "metadata.json"
    )
    metadata = {
        "model": cfg["model"],
        "model_version": cfg["model_version"],
        "experiment": cfg["experiment"],
        "run_date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": cfg,
        "baseline_summary": {
            "M1A_q_near": {k: (str(v) if k == "stopped_at" else v)
                           for k, v in trajectory_summary(m1a["rows"]).items()},
            "M0_q_near_T100": {"best_round": m0_best, "best_U": m0_u[m0_best],
                               "U_final": m0_u[-1]},
            "first_U_divergence_M1A_vs_M0": first_div,
            "alpha_first_state_divergence": state_div,
            "mc_boundary_flips_per_seed": dict(zip(seeds, boundary_flips)),
            "noise_floor": floor,
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m1a_adaptive_boundary/run_experiment.py",
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
