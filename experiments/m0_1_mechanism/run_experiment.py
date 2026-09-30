"""Runner for the M0.1 correction-strength time-scale experiment.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m0_1_mechanism/run_experiment.py [--config CONFIG.json]

Design (fixed in config.json):
  * deterministic trajectories for alpha in {0.0, 0.10, 0.25, 0.50, 0.75, 1.00}
    on the canonical M0-A configuration (q_near, a=0.50, K=100, T=30);
  * exact mechanism predictions from src/sand_m0/mechanism.py (closed-form
    state, near-mass recurrence, quadratic U(z), continuous optimum);
  * Monte Carlo spot checks (5 seeds, N=100000) for alpha in
    {0.25, 0.75, 1.00} to confirm the particle version does not
    systematically deviate from the deterministic mechanism.

This is NOT a parameter optimization: no alpha is declared "best".
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

from sand_m0.mechanism import (  # noqa: E402
    decay_factor,
    mechanism_constants,
    optimal_round,
    optimal_z,
    quadratic_coefficients,
    u_of_z,
)
from sand_m0.model import DISTRIBUTIONS, noise_floor_l2  # noqa: E402
from sand_m0.plotting import (  # noqa: E402
    plot_alpha_rounds,
    plot_mc_spotcheck,
    plot_rescaling,
)
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

    K, T, a, N = cfg["K"], cfg["T"], cfg["a"], cfg["N"]
    alphas = list(cfg["alphas_deterministic"])
    mc_alphas = list(cfg["alphas_mc_spotcheck"])
    seeds = list(cfg["mc_seeds"])
    started = time.time()

    print(f"M0.1 mechanism run — config: {config_path}")
    print(f"  a={a}  K={K}  T={T}  N={N}")
    print(f"  alphas (deterministic) = {alphas}")
    print(f"  alphas (MC spot check) = {mc_alphas}, seeds = {seeds}")

    dist = DISTRIBUTIONS["near"]
    C = mechanism_constants(dist, K, a)
    A_q, B_q, C_q = quadratic_coefficients(C)
    z_star = optimal_z(C)

    # shared theory curve U(z) on a dense z grid (identical for every alpha)
    z_grid = np.linspace(0.0, 1.0, 2001)
    u_disc = [A_q * z * z + B_q * z + C_q for z in z_grid]
    u_cont = [
        C.S_near * z * z
        + C.S_far * ((1.0 - 0.75 * z) / 0.25) ** 2
        - 1.0 / K
        for z in z_grid
    ]
    theory_csv = out_dir / "theory_u_of_z.csv"
    write_csv(
        theory_csv,
        ["z", "U_discrete_exact", "U_continuous"],
        [
            {"z": f"{z:.6f}", "U_discrete_exact": f"{ud:.12e}", "U_continuous": f"{uc:.12e}"}
            for z, ud, uc in zip(z_grid, u_disc, u_cont)
        ],
    )

    # ---------------- deterministic sweep ----------------
    det_curves: dict[float, dict] = {}
    traj_rows: list[dict] = []
    comparison_rows: list[dict] = []

    for alpha in alphas:
        rows = run_deterministic(dist, K=K, T=T, a=a, alpha=alpha)
        rho = decay_factor(alpha, C.Q_near)
        z_t = [rho**r["t"] for r in rows]
        u_sim = [r["U_L2"] for r in rows]
        u_theory = [u_of_z(z, C) for z in z_t]
        max_u_err = max(abs(us - ut) for us, ut in zip(u_sim, u_theory))
        max_m_err = max(abs(r["near_mass"] - C.Q_near * z) for r, z in zip(rows, z_t))

        for r, z, ut in zip(rows, z_t, u_theory):
            traj_rows.append(
                {
                    "alpha": alpha,
                    "t": r["t"],
                    "z_theory": f"{z:.12f}",
                    "near_mass": f"{r['near_mass']:.12f}",
                    "near_mass_theory": f"{C.Q_near * z:.12f}",
                    "U_L2": f"{r['U_L2']:.12e}",
                    "TV": f"{r['TV']:.12e}",
                    "U_theory": f"{ut:.12e}",
                }
            )

        if alpha == 0.0:
            best_t, best_u = 0, rows[0]["U_L2"]
            t_star = float("nan")
            note = "boundary: no operation, state never moves"
        else:
            best_t, best_u = best_round(rows)
            t_star = optimal_round(rho, z_star)
            note = ""

        comparison_rows.append(
            {
                "alpha": alpha,
                "rho_progression_factor": f"{rho:.6f}",
                "t_star_continuous": "" if t_star != t_star else f"{t_star:.4f}",
                "observed_best_round": best_t,
                "observed_best_U": f"{best_u:.12e}",
                "continuous_min_U_at_z_star": f"{u_of_z(z_star, C):.12e}",
                "discretization_gap_best_U_minus_U_zstar": f"{best_u - u_of_z(z_star, C):.12e}",
                "best_state_near_mass": f"{rows[best_t]['near_mass']:.10f}",
                "best_state_near_mass_theory_075_zstar": f"{C.Q_near * z_star:.10f}",
                "max_abs_near_mass_minus_theory": f"{max_m_err:.3e}",
                "max_abs_U_minus_quadratic": f"{max_u_err:.3e}",
                "note": note,
            }
        )
        det_curves[alpha] = {"u": u_sim, "z": z_t, "best_t": best_t}

    alpha_csv = out_dir / "alpha_comparison.csv"
    write_csv(
        alpha_csv,
        ["alpha", "rho_progression_factor", "t_star_continuous", "observed_best_round",
         "observed_best_U", "continuous_min_U_at_z_star",
         "discretization_gap_best_U_minus_U_zstar", "best_state_near_mass",
         "best_state_near_mass_theory_075_zstar", "max_abs_near_mass_minus_theory",
         "max_abs_U_minus_quadratic", "note"],
        comparison_rows,
    )

    traj_csv = out_dir / "deterministic_trajectories.csv"
    write_csv(
        traj_csv,
        ["alpha", "t", "z_theory", "near_mass", "near_mass_theory", "U_L2", "TV",
         "U_theory"],
        traj_rows,
    )

    # ---------------- MC spot checks ----------------
    floor = noise_floor_l2(N, K)
    mc_rows: list[dict] = []
    mc_max_diff = {}
    mc_mean_max_diff = {}
    for alpha in mc_alphas:
        rho = decay_factor(alpha, C.Q_near)
        per_t: dict[int, list[float]] = {}
        for seed in seeds:
            run = run_mc(dist, K=K, T=T, a=a, alpha=alpha, N=N, seed=seed)
            for r in run:
                z = rho**r["t"]
                ut = u_of_z(z, C)
                per_t.setdefault(r["t"], []).append(r["U_L2"])
                mc_rows.append(
                    {
                        "alpha": alpha,
                        "seed": seed,
                        "t": r["t"],
                        "z": f"{z:.12f}",
                        "near_mass": f"{r['near_mass']:.10f}",
                        "U_L2": f"{r['U_L2']:.12e}",
                        "U_theory": f"{ut:.12e}",
                        "U_minus_theory": f"{r['U_L2'] - ut:.12e}",
                    }
                )
        # per-seed scatter (expected: multinomial sampling noise)
        mc_max_diff[alpha] = max(
            abs(float(row["U_minus_theory"])) for row in mc_rows if row["alpha"] == alpha
        )
        # systematic-deviation check: MC mean per round vs theory + noise floor.
        # E[U_mc] = U(z_t) + (1 - 1/K)/N for a multinomial histogram.
        mc_mean_max_diff[alpha] = max(
            abs(np.mean(per_t[t]) - (u_of_z(rho**t, C) + floor))
            for t in sorted(per_t)
        )

    mc_csv = out_dir / "mc_spotcheck.csv"
    write_csv(
        mc_csv,
        ["alpha", "seed", "t", "z", "near_mass", "U_L2", "U_theory", "U_minus_theory"],
        mc_rows,
    )

    # ---------------- figures ----------------
    fig1 = out_dir / "fig_alpha_uniformity_vs_round.png"
    plot_alpha_rounds(det_curves, path=fig1, T=T)
    fig2 = out_dir / "fig_rescaling_u_vs_z.png"
    plot_rescaling(
        det_curves, z_grid=z_grid, u_disc=u_disc, u_cont=u_cont,
        z_star=z_star, constants=C, path=fig2,
    )
    fig3 = out_dir / "fig_mc_spotcheck.png"
    plot_mc_spotcheck(
        mc_csv=mc_csv, z_grid=z_grid, u_disc=u_disc, noise_floor=floor, path=fig3
    )

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
        "theory_summary": {
            "Q_near": C.Q_near,
            "Q_far": C.Q_far,
            "S_near": C.S_near,
            "S_far": C.S_far,
            "baseline_rho_15_over_16": decay_factor(0.25, C.Q_near),
            "z_star_discrete": z_star,
            "z_star_correction_from_075": z_star - 0.75,
            "t_star_baseline": optimal_round(decay_factor(0.25, C.Q_near), z_star),
            "U_quadratic_A_B_C": [A_q, B_q, C_q],
            "continuous_min_U_1_over_600": 1.0 / 600.0,
        },
        "mc_spotcheck_max_abs_U_minus_theory": {
            str(k): f"{v:.3e}" for k, v in mc_max_diff.items()
        },
        "mc_spotcheck_max_abs_mean_minus_theory_plus_floor": {
            str(k): f"{v:.3e}" for k, v in mc_mean_max_diff.items()
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m0_1_mechanism/run_experiment.py",
        "outputs": {name: sha256_of(out_dir / name) for name in output_files},
    }
    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\n=== alpha comparison (deterministic) ===")
    for row in comparison_rows:
        print(
            f"  alpha={row['alpha']:<5} rho={row['rho_progression_factor']}  "
            f"t*={row['t_star_continuous'] or 'n/a':>7}  "
            f"best_round={row['observed_best_round']:>2}  "
            f"best_U={float(row['observed_best_U']):.6e}  "
            f"gap_to_U(z*)={float(row['discretization_gap_best_U_minus_U_zstar']):.2e}"
        )
    print("\n=== MC spot check per alpha ===")
    for k in mc_max_diff:
        print(f"  alpha={k}: max |U_mc - U_theory| (single seed) = {mc_max_diff[k]:.3e}; "
              f"max |MC mean - (theory + noise floor)| = {mc_mean_max_diff[k]:.3e}")
    print(f"\noutputs written to: {out_dir} ({len(output_files) + 1} files)")
    print(f"runtime: {metadata['runtime_seconds']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
