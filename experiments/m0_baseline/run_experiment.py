"""Runner for the three official M0 experiments (EXP-M0-A/B/C).

One command regenerates every CSV, figure and the metadata file:

    python experiments/m0_baseline/run_experiment.py [--config CONFIG.json]

Fixed parameters, seeds and the metric definitions live in config.json next
to this script. Outputs are written to <config dir>/results/. No parameter
tuning happens here; changing config.json is a new experiment and must be
declared as such in the research log.
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

from sand_m0.model import DISTRIBUTIONS, noise_floor_l2  # noqa: E402
from sand_m0.plotting import plot_det_vs_mc, plot_experiment  # noqa: E402
from sand_m0.simulate import best_round, run_deterministic, run_mc  # noqa: E402


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

    commit = run(["rev-parse", "HEAD"])
    dirty = bool(run(["status", "--porcelain"]))
    return {"commit": commit, "working_tree_dirty_at_run_time": dirty}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path(__file__).resolve().parent / "config.json",
        help="path to the experiment config JSON (default: config.json next to this script)",
    )
    args = parser.parse_args(argv)

    config_path = args.config.resolve()
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    out_dir = config_path.parent / cfg["output_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)

    K, T, a, alpha, N = cfg["K"], cfg["T"], cfg["a"], cfg["alpha"], cfg["N"]
    seeds = list(cfg["seeds"])
    started = time.time()

    print(f"M0 baseline run — config: {config_path}")
    print(f"  a={a}  alpha={alpha}  K={K}  T={T}  N={N}  seeds={len(seeds)}")

    det_curves: dict[str, dict] = {}
    summary_for_report: list[str] = []

    for exp in cfg["experiments"]:
        exp_id = exp["id"]
        dist = DISTRIBUTIONS[exp["distribution"]]
        print(f"\n[{exp_id}] distribution={dist.name}: {exp['question']}")

        det_rows = run_deterministic(dist, K=K, T=T, a=a, alpha=alpha)
        det_best, det_best_u = best_round(det_rows)

        mc_per_seed: list[list[dict]] = []
        for seed in seeds:
            mc_per_seed.append(run_mc(dist, K=K, T=T, a=a, alpha=alpha, N=N, seed=seed))

        ts = [row["t"] for row in det_rows]
        mc_u = np.array([[row["U_L2"] for row in run] for run in mc_per_seed])
        mc_tv = np.array([[row["TV"] for row in run] for run in mc_per_seed])
        mc_mean, mc_std = mc_u.mean(axis=0), mc_u.std(axis=0, ddof=1)
        mc_tv_mean = mc_tv.mean(axis=0)
        mc_best, mc_best_u = min(zip(ts, mc_mean), key=lambda pair: pair[1])

        floor = noise_floor_l2(N, K)

        det_csv = out_dir / f"{exp_id}_deterministic.csv"
        write_csv(
            det_csv,
            ["t", "near_mass", "removed_mass", "U_L2", "TV"],
            [
                {
                    "t": row["t"],
                    "near_mass": f"{row['near_mass']:.10f}",
                    "removed_mass": f"{row['removed_mass']:.10f}",
                    "U_L2": f"{row['U_L2']:.12e}",
                    "TV": f"{row['TV']:.12e}",
                }
                for row in det_rows
            ],
        )

        seed_csv = out_dir / f"{exp_id}_mc_per_seed.csv"
        write_csv(
            seed_csv,
            ["seed", "t", "near_mass", "removed_frac", "U_L2", "TV"],
            [
                {
                    "seed": row["seed"],
                    "t": row["t"],
                    "near_mass": f"{row['near_mass']:.10f}",
                    "removed_frac": f"{row['removed_frac']:.10f}",
                    "U_L2": f"{row['U_L2']:.12e}",
                    "TV": f"{row['TV']:.12e}",
                }
                for run in mc_per_seed
                for row in run
            ],
        )

        summary_csv = out_dir / f"{exp_id}_mc_summary.csv"
        write_csv(
            summary_csv,
            ["t", "mc_mean_U_L2", "mc_std_U_L2", "mc_mean_TV", "det_U_L2",
             "mc_mean_minus_det_U_L2", "noise_floor_U_L2"],
            [
                {
                    "t": t,
                    "mc_mean_U_L2": f"{mc_mean[i]:.12e}",
                    "mc_std_U_L2": f"{mc_std[i]:.12e}",
                    "mc_mean_TV": f"{mc_tv_mean[i]:.12e}",
                    "det_U_L2": f"{det_rows[i]['U_L2']:.12e}",
                    "mc_mean_minus_det_U_L2": f"{mc_mean[i] - det_rows[i]['U_L2']:.12e}",
                    "noise_floor_U_L2": f"{floor:.12e}",
                }
                for i, t in enumerate(ts)
            ],
        )

        best_csv = out_dir / f"{exp_id}_mc_seed_best.csv"
        write_csv(
            best_csv,
            ["seed", "t_best", "U_best", "U_at_t0", "U_at_T"],
            [
                {
                    "seed": run[0]["seed"],
                    "t_best": t_best,
                    "U_best": f"{u_best:.12e}",
                    "U_at_t0": f"{run[0]['U_L2']:.12e}",
                    "U_at_T": f"{run[-1]['U_L2']:.12e}",
                }
                for run in mc_per_seed
                for t_best, u_best in [best_round(run)]
            ],
        )

        fig_path = out_dir / f"fig_{exp_id}_uniformity.png"
        plot_experiment(
            det_rows, mc_mean, mc_std,
            exp_id=exp_id, dist_label=dist.label, noise_floor=floor,
            det_best=det_best, mc_best=int(mc_best), path=fig_path,
        )

        det_curves[exp_id] = {
            "label": dist.label,
            "det_u": [row["U_L2"] for row in det_rows],
            "mc_mean": mc_mean,
        }

        from collections import Counter

        best_counts = Counter(best_round(run)[0] for run in mc_per_seed)
        summary_for_report.append(
            f"{exp_id} [{dist.name}]:\n"
            f"  det:  U(0)={det_rows[0]['U_L2']:.6e}  t_best={det_best} "
            f"(U={det_best_u:.6e})  U(T)={det_rows[-1]['U_L2']:.6e}\n"
            f"  det:  near_mass: t0={det_rows[0]['near_mass']:.4f} "
            f"t5={det_rows[5]['near_mass']:.4f} t10={det_rows[10]['near_mass']:.4f} "
            f"t30={det_rows[-1]['near_mass']:.4f}\n"
            f"  mc:   mean best t={int(mc_best)} (U={mc_best_u:.6e})  "
            f"U(0) mean={mc_mean[0]:.6e}  U(T) mean={mc_mean[-1]:.6e}\n"
            f"  mc:   per-seed t_best distribution: {dict(sorted(best_counts.items()))}\n"
            f"  max |mc mean - det| U: {np.max(np.abs(mc_mean - np.array(det_curves[exp_id]['det_u']))):.3e}"
        )

    combined_fig = out_dir / "fig_det_vs_mc.png"
    plot_det_vs_mc(det_curves, path=combined_fig, N=N, n_seeds=len(seeds))

    output_files = sorted(
        p.name for p in out_dir.iterdir() if p.is_file() and p.name != "metadata.json"
    )
    metadata = {
        "model": cfg["model"],
        "model_version": cfg["model_version"],
        "run_date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": cfg,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m0_baseline/run_experiment.py",
        "outputs": {name: sha256_of(out_dir / name) for name in output_files},
    }
    metadata_path = out_dir / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\n=== summary ===")
    for block in summary_for_report:
        print(block)
    print(f"\noutputs written to: {out_dir} ({len(output_files)} files)")
    print(f"runtime: {metadata['runtime_seconds']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
