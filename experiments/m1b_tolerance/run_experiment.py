"""Runner for M1B: perceptual tolerance / deadband stopping.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m1b_tolerance/run_experiment.py

Only the stopping rule changes relative to M1A: stop when D_max <= epsilon
(finite tolerance) instead of D_max <= 1e-12. Because M1B follows the exact
M1A reference trajectory until it stops, the stopping time is a FIRST
PASSAGE of the reference D_max sequence through epsilon — so the reference
run (epsilon = 0, implemented with the canonical 1e-12 tolerance) predicts
every finite-epsilon stopping time before any M1B run. The runner verifies
this prediction against actual canonical runs of each (K, epsilon) case.

Frozen design: epsilons {0, 0.005, 0.01, 0.02, 0.04} x K {50, 100, 200, 400},
T_max = 5000, alpha = 0.25, q_near. Not a parameter optimization.
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
from sand_m0.model import DISTRIBUTIONS, noise_floor_l2  # noqa: E402
from sand_m0.plotting import (  # noqa: E402
    plot_m1b_reference,
    plot_stopping_vs_epsilon,
    plot_regret_vs_epsilon,
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
    K_list = list(cfg["K_list"])
    eps_entries = list(cfg["epsilons"])
    dense_grid = np.geomspace(1e-12, 0.06, 200)
    started = time.time()

    print(f"M1B run — config: {config_path}")
    print(f"  epsilons = {[e['epsilon'] for e in eps_entries]}, K = {K_list}, "
          f"T_max = {T_max}, alpha = {alpha}")

    # ---------------- 1. exact reference runs (epsilon = 0) per K ----------------
    references: dict[int, dict] = {}
    for K in K_list:
        t0 = time.time()
        references[K] = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K,
                                              T_max=T_max, alpha=alpha,
                                              stop_tolerance=1e-12)
        Ud = np.array([K * r["U_L2"] for r in references[K]["rows"]])
        best_t = int(Ud.argmin())
        references[K]["U_density"] = Ud
        references[K]["U_raw"] = np.array([r["U_L2"] for r in references[K]["rows"]])
        references[K]["D_max_seq"] = np.array([r["D_max"] for r in references[K]["rows"]])
        references[K]["a_seq"] = np.array([r["a_t"] for r in references[K]["rows"]])
        references[K]["near_seq"] = np.array([r["near_mass_05"] for r in references[K]["rows"]])
        references[K]["t_U_best"] = best_t
        references[K]["U_density_best"] = float(Ud[best_t])
        print(f"  [reference K={K}] stopped_at={references[K]['stopped_at']} "
              f"({time.time() - t0:.1f}s); t_U_best={best_t} "
              f"U_density_best={Ud[best_t]:.6e}")

    # dense first-passage table per K (monotonicity verification, no re-simulation)
    mono_rows: list[dict] = []
    for K in K_list:
        Dm = references[K]["D_max_seq"]
        ts = []
        for eps in dense_grid:
            hit = np.nonzero(Dm <= eps)[0]
            ts.append(int(hit[0]) if len(hit) else -1)
        non_increasing = all(
            ts[i] == -1 or ts[i + 1] == -1 or ts[i] >= ts[i + 1] for i in range(len(ts) - 1)
        )
        references[K]["dense_first_passage"] = ts
        mono_rows.append({"K": K, "dense_grid_non_increasing": non_increasing,
                          "n_grid_points": len(dense_grid)})
        print(f"  [monotonicity K={K}] t_stop non-increasing over 200-pt eps grid: {non_increasing}")
    mono_csv = out_dir / "first_passage_monotonicity.csv"
    write_csv(mono_csv, list(mono_rows[0]), mono_rows)

    # ---------------- 2. epsilon x K canonical runs ----------------
    case_rows: list[dict] = []
    for K in K_list:
        ref = references[K]
        for entry in eps_entries:
            eps = entry["epsilon"]
            tol = entry["implemented_tolerance"]
            run = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=T_max,
                                        alpha=alpha, stop_tolerance=tol)
            stopped = run["stopped_at"]
            rows = run["rows"]
            # first-passage prediction from the reference sequence
            hit = np.nonzero(ref["D_max_seq"] <= tol)[0]
            predicted = int(hit[0]) if len(hit) else -1
            if stopped is not None:
                st_row = rows[stopped]
                u_stop = rows[stopped]["U_L2"]
                ud_stop = K * u_stop
                a_before = rows[stopped - 1]["a_t"] if stopped >= 1 else float("nan")
                case_rows.append({
                    "K": K,
                    "epsilon": eps,
                    "implemented_tolerance": tol,
                    "predicted_t_stop_from_reference": predicted,
                    "matches_reference_prediction": (stopped == predicted),
                    "stopping_round": stopped,
                    "D_max_at_stop": f"{st_row['D_max']:.12e}",
                    "U_density_at_stop": f"{ud_stop:.12e}",
                    "U_raw_at_stop": f"{u_stop:.12e}",
                    "near_half_mass_at_stop": f"{st_row['near_mass_05']:.10f}",
                    "boundary_before_stop": f"{a_before:.4f}" if stopped >= 1 else "",
                    "actions_executed": stopped,
                    "delta_t_vs_U_best": stopped - ref["t_U_best"],
                    "uniformity_regret": f"{(ud_stop - ref['U_density_best']) / ref['U_density_best']:.6f}",
                })
            else:
                case_rows.append({
                    "K": K, "epsilon": eps, "implemented_tolerance": tol,
                    "predicted_t_stop_from_reference": predicted,
                    "matches_reference_prediction": (predicted == -1),
                    "stopping_round": "none",
                    "D_max_at_stop": "", "U_density_at_stop": "", "U_raw_at_stop": "",
                    "near_half_mass_at_stop": "", "boundary_before_stop": "",
                    "actions_executed": f">{T_max}",
                    "delta_t_vs_U_best": "", "uniformity_regret": "",
                })
    cases_csv = out_dir / "epsilon_k_summary.csv"
    write_csv(cases_csv, list(case_rows[0]), case_rows)

    print("\n=== epsilon x K results ===")
    print("  K    eps    t_stop  pred  match  U_dens@stop  regret   near")
    for r in case_rows:
        if r["stopping_round"] != "none":
            print(f"  {r['K']:>3} {r['epsilon']:>6} {r['stopping_round']:>6} "
                  f"{r['predicted_t_stop_from_reference']:>5} "
                  f"{str(r['matches_reference_prediction']):>5}  "
                  f"{float(r['U_density_at_stop']):.5e}  "
                  f"{float(r['uniformity_regret']):>7.3f}  "
                  f"{float(r['near_half_mass_at_stop']):.4f}")
        else:
            print(f"  {r['K']:>3} {r['epsilon']:>6}  no stop within T_max")

    # ---------------- 3. controls (sanity) ----------------
    control_rows = []
    for name in cfg["controls"]["cases"]:
        run = run_m1a_deterministic(DISTRIBUTIONS[name], K=100, T_max=5,
                                    alpha=alpha, stop_tolerance=0.01)
        r0 = run["rows"][0]
        control_rows.append({
            "distribution": f"q_{name}",
            "K": 100,
            "epsilon": 0.01,
            "D_max_t0": f"{r0['D_max']:.6e}",
            "stopped_at": str(run["stopped_at"]),
            "sweeps_executed": sum(1 for r in run["rows"] if r["active"]),
        })
        print(f"[control q_{name}] stopped_at={run['stopped_at']} "
              f"D_max(0)={r0['D_max']:.2e}")
    ctrl_csv = out_dir / "controls.csv"
    write_csv(ctrl_csv, list(control_rows[0]), control_rows)

    # ---------------- 4. figures ----------------
    fig1 = out_dir / "fig_stopping_vs_epsilon.png"
    plot_stopping_vs_epsilon(case_rows, K_list, eps_entries, path=fig1)
    fig2 = out_dir / "fig_regret_vs_epsilon.png"
    plot_regret_vs_epsilon(case_rows, K_list, eps_entries, path=fig2)
    fig3 = out_dir / "fig_m1b_reference_D_max.png"
    plot_m1b_reference(references, K_list, eps_entries, path=fig3)

    # ---------------- 5. metadata ----------------
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
            "reference_stopped_at_per_K": {str(K): str(references[K]["stopped_at"]) for K in K_list},
            "t_U_best_per_K": {str(K): references[K]["t_U_best"] for K in K_list},
            "all_cases_match_reference_first_passage": all(
                r["matches_reference_prediction"] for r in case_rows
            ),
            "monotonicity_verified_per_K": {str(r["K"]): r["dense_grid_non_increasing"]
                                            for r in mono_rows},
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m1b_tolerance/run_experiment.py",
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
