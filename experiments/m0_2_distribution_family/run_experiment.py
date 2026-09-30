"""Runner for M0.2 Phase 2A: throw distribution family.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m0_2_distribution_family/run_experiment.py

Order of operations inside this runner encodes the epistemics of the phase:
  1. THEORY FIRST — for every frozen family member, compute the
     pre-simulation prediction table (Q, shape concentrations, predicted
     existence of an interior optimum, predicted z*, predicted continuous
     t*, predicted discrete best round) and write predictions.csv;
  2. only then run the canonical deterministic simulator per distribution;
  3. compare prediction vs observation (comparison.csv);
  4. alpha time-rescaling spot checks for a few representative q
     (deterministic; exactness against the shared quadratic);
  5. Monte Carlo spot checks (5 seeds) to confirm the particle version
     fluctuates around the deterministic prediction.

No alpha grid, no parameter optimization: baseline alpha = 0.25 for the
family comparison; alpha spot checks use the pre-declared set {0.1, 0.5, 1.0}
on three representative distributions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from sand_m0.mechanism import (  # noqa: E402
    analytic_trajectory,
    decay_factor,
    mechanism_constants,
    optimal_round,
    optimal_z,
    quadratic_coefficients,
    u_of_z,
)
from sand_m0.model import noise_floor_l2, power_distribution  # noqa: E402
from sand_m0.plotting import (  # noqa: E402
    plot_family_shapes,
    plot_prediction_vs_observed,
    plot_u_vs_round_family,
    plot_u_vs_z_family,
)
from sand_m0.simulate import run_deterministic, run_mc  # noqa: E402


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


def predict(p: float, *, K: int, a: float, alpha: float, T: int) -> dict:
    """Pre-simulation theoretical prediction for one family member.

    Uses only the model definition: exact discrete bin probabilities of q_p,
    the exact zone constants, the closed-form path z(t), the exact quadratic
    U(z) and its minimizer. No simulation is involved.
    """
    dist = power_distribution(p)
    C = mechanism_constants(dist, K, a)
    rho = decay_factor(alpha, C.Q_near)
    z_star = optimal_z(C)
    c_n = C.S_near / C.Q_near
    c_f = C.S_far / C.Q_far
    u_cont = 2.0 ** (-(p + 1.0))
    z_star_cont = 2.0 * (1.0 - u_cont) / (3.0 - 4.0 * u_cont)
    t_star_cont = (
        optimal_round(rho, z_star_cont) if z_star_cont < 1.0 else float("nan")
    )
    A, B, Cc = quadratic_coefficients(C)
    u_path = [A * (rho**t) ** 2 + B * (rho**t) + Cc for t in range(T + 1)]
    best_t = min(range(T + 1), key=lambda t: u_path[t])
    margins = sorted(abs(u_path[t] - min(u_path)) for t in range(T + 1))
    interior = z_star < 1.0 - 1e-9  # tolerance: the exact boundary (z*=1) is not interior
    return {
        "p": p,
        "Q_near": C.Q_near,
        "S_near": C.S_near,
        "S_far": C.S_far,
        "c_n": c_n,
        "c_f": c_f,
        "concentration_ratio_c_n_over_c_f": c_n / c_f,
        "predicted_interior_optimum": interior,
        "z_star_discrete": z_star,
        "z_star_continuous_closed_form": z_star_cont,
        "rho_progression_factor": rho,
        "t_star_continuous": t_star_cont,
        "predicted_best_round_discrete": best_t,
        "predicted_best_U": u_path[best_t],
        "predicted_U_at_t0": u_path[0],
        "predicted_U_worse_at_end_than_start": u_path[-1] > u_path[0],
        "best_U_margin_to_runner_up": margins[1] - margins[0],
        "_dist": dist,
        "_constants": C,
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

    K, T, a, N = cfg["K"], cfg["T"], cfg["a"], cfg["N"]
    alpha0 = cfg["baseline_alpha"]
    ps = list(cfg["family"]["cases_p"])
    seeds = list(cfg["mc_spotcheck"]["seeds"])
    started = time.time()

    print(f"M0.2 Phase 2A run — config: {config_path}")
    print(f"  family: {cfg['family']['definition']}")
    print(f"  cases p = {ps}, baseline alpha = {alpha0}, K = {K}, T = {T}")

    # ---------------- 1. theory first: pre-simulation predictions ----------------
    preds = [predict(p, K=K, a=a, alpha=alpha0, T=T) for p in ps]
    pred_csv = out_dir / "predictions.csv"
    pred_fields = [k for k in preds[0] if not k.startswith("_")]
    write_csv(pred_csv, pred_fields, [{k: row[k] for k in pred_fields} for row in preds])
    print("\n=== predictions written BEFORE simulation (predictions.csv) ===")
    for row in preds:
        print(
            f"  p={row['p']:>5}: interior={str(row['predicted_interior_optimum']):>5}  "
            f"z*={row['z_star_discrete']:.5f}  t*={row['t_star_continuous']:>7.3f}  "
            f"pred_best_round={row['predicted_best_round_discrete']:>2}"
        )

    # ---------------- 2. canonical simulation per distribution ----------------
    traj_rows: list[dict] = []
    comp_rows: list[dict] = []
    curves: dict[float, dict] = {}

    for pred in preds:
        p = pred["p"]
        dist = pred["_dist"]
        C = pred["_constants"]
        rows = run_deterministic(dist, K=K, T=T, a=a, alpha=alpha0)
        rho = decay_factor(alpha0, C.Q_near)

        max_u_err = max(abs(r["U_L2"] - u_of_z(rho**r["t"], C)) for r in rows)
        best_t = min(range(T + 1), key=lambda t: rows[t]["U_L2"])
        best_u = rows[best_t]["U_L2"]
        z_best_obs = rows[best_t]["near_mass"] / C.Q_near

        for r in rows:
            traj_rows.append(
                {
                    "p": p,
                    "t": r["t"],
                    "z": f"{rho**r['t']:.12f}",
                    "near_mass": f"{r['near_mass']:.12f}",
                    "U_L2": f"{r['U_L2']:.12e}",
                    "TV": f"{r['TV']:.12e}",
                    "U_theory": f"{u_of_z(rho**r['t'], C):.12e}",
                }
            )
        comp_rows.append(
            {
                "p": p,
                "role": cfg["family"]["roles"].get(str(p), ""),
                "predicted_interior_optimum": pred["predicted_interior_optimum"],
                "observed_improves_then_worsens": best_t > 0 and rows[-1]["U_L2"] > best_u,
                "predicted_best_round": pred["predicted_best_round_discrete"],
                "observed_best_round": best_t,
                "predicted_z_star": f"{pred['z_star_discrete']:.6f}",
                "observed_best_state_z": f"{z_best_obs:.6f}",
                "predicted_best_U": f"{pred['predicted_best_U']:.12e}",
                "observed_best_U": f"{best_u:.12e}",
                "observed_minus_predicted_best_U": f"{best_u - pred['predicted_best_U']:.3e}",
                "max_abs_U_minus_quadratic": f"{max_u_err:.3e}",
                "observed_U_t0": f"{rows[0]['U_L2']:.12e}",
                "observed_U_T": f"{rows[-1]['U_L2']:.12e}",
            }
        )
        curves[p] = {"u": [r["U_L2"] for r in rows], "z": [rho**r["t"] for r in rows],
                     "best_t": best_t}

    traj_csv = out_dir / "deterministic_trajectories.csv"
    write_csv(traj_csv, ["p", "t", "z", "near_mass", "U_L2", "TV", "U_theory"], traj_rows)
    comp_csv = out_dir / "comparison.csv"
    write_csv(comp_csv, list(comp_rows[0]), comp_rows)

    print("\n=== prediction vs observation ===")
    for row in comp_rows:
        print(
            f"  p={row['p']:>5}: pred_round={row['predicted_best_round']:>2}  "
            f"obs_round={row['observed_best_round']:>2}  "
            f"interior pred/obs = {row['predicted_interior_optimum']}/"
            f"{row['observed_improves_then_worsens']}  "
            f"z* pred/obs = {float(row['predicted_z_star']):.5f}/"
            f"{float(row['observed_best_state_z']):.5f}  "
            f"max|U-U(z)|={row['max_abs_U_minus_quadratic']}"
        )

    # ---------------- 3. alpha time-rescaling spot checks ----------------
    alpha_rows: list[dict] = []
    alpha_traj: dict[float, dict[float, tuple[list[float], list[float]]]] = {}
    for p in cfg["alpha_spotcheck"]["cases_p"]:
        pred = next(r for r in preds if r["p"] == p)
        C = pred["_constants"]
        for alpha in cfg["alpha_spotcheck"]["alphas"]:
            dist = power_distribution(p)
            rows = run_deterministic(dist, K=K, T=T, a=a, alpha=alpha)
            rho = decay_factor(alpha, C.Q_near)
            max_u_err = max(abs(r["U_L2"] - u_of_z(rho**r["t"], C)) for r in rows)
            best_t = min(range(T + 1), key=lambda t: rows[t]["U_L2"])
            best_u = rows[best_t]["U_L2"]
            z_star = optimal_z(C)
            t_star = optimal_round(rho, z_star) if z_star < 1 else float("nan")
            alpha_traj.setdefault(p, {})[alpha] = (
                [rho**r["t"] for r in rows],
                [r["U_L2"] for r in rows],
            )
            alpha_rows.append(
                {
                    "p": p,
                    "alpha": alpha,
                    "rho_progression_factor": f"{rho:.6f}",
                    "z_star_discrete": f"{z_star:.6f}",
                    "t_star_continuous": "" if t_star != t_star else f"{t_star:.4f}",
                    "observed_best_round": best_t,
                    "observed_best_U": f"{best_u:.12e}",
                    "continuous_min_U_at_z_star": f"{u_of_z(z_star, C):.12e}",
                    "max_abs_U_minus_quadratic": f"{max_u_err:.3e}",
                    "note": "time-rescaling spot check" if z_star < 1 else
                            "boundary/counterexample: t=0 optimal for every alpha",
                }
            )
    alpha_csv = out_dir / "alpha_rescaling.csv"
    write_csv(alpha_csv, list(alpha_rows[0]), alpha_rows)

    # ---------------- 4. Monte Carlo spot checks ----------------
    floor = noise_floor_l2(N, K)
    mc_rows: list[dict] = []
    mc_mean_dev = {}
    for p in cfg["mc_spotcheck"]["cases_p"]:
        pred = next(r for r in preds if r["p"] == p)
        C = pred["_constants"]
        rho = decay_factor(alpha0, C.Q_near)
        per_t: dict[int, list[float]] = {}
        for seed in seeds:
            run = run_mc(power_distribution(p), K=K, T=T, a=a, alpha=alpha0,
                         N=N, seed=seed)
            for r in run:
                z = rho**r["t"]
                per_t.setdefault(r["t"], []).append(r["U_L2"])
                mc_rows.append(
                    {
                        "p": p,
                        "seed": seed,
                        "t": r["t"],
                        "z": f"{z:.12f}",
                        "U_L2": f"{r['U_L2']:.12e}",
                        "U_theory": f"{u_of_z(z, C):.12e}",
                        "U_minus_theory": f"{r['U_L2'] - u_of_z(z, C):.12e}",
                    }
                )
        mc_mean_dev[p] = max(
            abs(np.mean(per_t[t]) - (u_of_z(rho**t, C) + floor)) for t in sorted(per_t)
        )
    mc_csv = out_dir / "mc_spotcheck.csv"
    write_csv(mc_csv, ["p", "seed", "t", "z", "U_L2", "U_theory", "U_minus_theory"],
              mc_rows)

    # ---------------- 5. figures ----------------
    shapes_fig = out_dir / "fig_family_shapes.png"
    plot_family_shapes(
        {p: power_distribution(p) for p in ps}, a=a, path=shapes_fig
    )
    rounds_fig = out_dir / "fig_uniformity_vs_round.png"
    plot_u_vs_round_family(curves, alpha=alpha0, path=rounds_fig)
    z_fig = out_dir / "fig_u_vs_z_family.png"
    plot_u_vs_z_family(
        curves,
        theory={p: quadratic_coefficients(pred["_constants"]) for p, pred in zip(ps, preds)},
        alpha_trajectories=alpha_traj,
        path=z_fig,
    )
    pred_fig = out_dir / "fig_prediction_vs_observed.png"
    plot_prediction_vs_observed(comp_rows, alpha_rows, path=pred_fig)

    # ---------------- 6. metadata ----------------
    output_files = sorted(
        pth.name for pth in out_dir.iterdir() if pth.is_file() and pth.name != "metadata.json"
    )
    metadata = {
        "model": cfg["model"],
        "model_version": cfg["model_version"],
        "experiment": cfg["experiment"],
        "family": cfg["family"],
        "run_date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": cfg,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "mc_spotcheck_max_abs_mean_minus_theory_plus_floor": {
            str(p): f"{v:.3e}" for p, v in mc_mean_dev.items()
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m0_2_distribution_family/run_experiment.py",
        "outputs": {name: sha256_of(out_dir / name) for name in output_files},
    }
    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\n=== MC spot check per p (baseline alpha) ===")
    for p, v in mc_mean_dev.items():
        print(f"  p={p}: max |MC mean - (U(z)+noise floor)| = {v:.3e} (floor {floor:.2e})")
    print(f"\noutputs written to: {out_dir} ({len(output_files) + 1} files)")
    print(f"runtime: {metadata['runtime_seconds']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
