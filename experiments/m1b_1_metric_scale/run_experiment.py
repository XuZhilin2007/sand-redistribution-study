"""Runner for M1B.1: uniformity metric & observation-scale audit.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m1b_1_metric_scale/run_experiment.py

Diagnostics only — the M1A/M1B controller, stopping rule and dynamics are
untouched. The exact M1A reference trajectories (identical canonical
dynamics, regenerated for state histories) are evaluated at FIXED coarse
observation scales B in {5,10,25,50} (every B divides every K, so the
aggregation P_b = sum of K/B consecutive bins is exact) alongside the
fine-grid metric U_K = K*sum(p_i - 1/K)^2.

For every (K, B): the reference optimum (t_best, U_best). For every
(K, epsilon): the M1B stopping state evaluated at every B — absolute gap,
residual regret, improvement loss / retention. Regression anchors: the
stopping times must reproduce the committed M1B results, and the K=100
trajectory must reproduce the committed M1A.1 trajectory.
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
    plot_t_best_vs_K,
    plot_retention_vs_epsilon,
    plot_regret_vs_loss,
    plot_coarse_profiles,
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
    tol = cfg["stop_tolerance_exact"]
    K_list = list(cfg["K_list"])
    B_list = list(cfg["B_list"])
    epsilons = list(cfg["epsilons"])
    started = time.time()

    print(f"M1B.1 metric-scale audit — config: {config_path}")
    print(f"  K={K_list}, B={B_list}, epsilons={epsilons}, T_max={T_max}")

    # ---------------- 1. exact reference runs + coarse evaluation ----------------
    runs: dict[int, dict] = {}
    U_by_B: dict[int, dict[int, np.ndarray]] = {}
    for K in K_list:
        t0 = time.time()
        run = run_m1a_deterministic(DISTRIBUTIONS["near"], K=K, T_max=T_max,
                                    alpha=alpha, stop_tolerance=tol)
        states = np.asarray(run["states"], dtype=float)  # (T+1, K)
        run["states_arr"] = states
        U_by_B[K] = {}
        for B in B_list + [K]:  # B=K reduces to the fine-grid U_density exactly
            P = states.reshape(states.shape[0], B, K // B).sum(axis=2)
            U_by_B[K][B] = B * ((P - 1.0 / B) ** 2).sum(axis=1)
        runs[K] = run
        print(f"  [K={K}] reference done ({time.time() - t0:.1f}s), "
              f"stopped_at={run['stopped_at']}")

    # ---------------- 2. fixed-scale reference optima ----------------
    opt_rows: list[dict] = []
    for K in K_list:
        for B in sorted(set(B_list) | {K}):  # B=K is the fine-grid metric; dedupe
            u = U_by_B[K][B]
            t_best = int(u.argmin())
            opt_rows.append({
                "K": K,
                "B": B,
                "scale_label": "fine (B=K)" if B == K else f"coarse B={B}",
                "t_best": t_best,
                "U_best": f"{u[t_best]:.12e}",
                "U_initial": f"{u[0]:.12e}",
                "improvement_available": f"{u[0] - u[t_best]:.12e}",
            })
    opt_csv = out_dir / "fixed_scale_optimum.csv"
    write_csv(opt_csv, list(opt_rows[0]), opt_rows)
    print("\n=== reference optima (t_best, U_best) per (K, B) ===")
    for r in opt_rows:
        print(f"  K={r['K']:>3} {r['scale_label']:<12} t_best={r['t_best']:>4} "
              f"U_best={float(r['U_best']):.5e}")

    # ---------------- 3. stop-quality across (K, epsilon, B) ----------------
    # stopping states are read off the reference trajectory (M1B stops are
    # first passages of the same sequence — verified against M1B CSV below)
    stop_rows: list[dict] = []
    for K in K_list:
        ref = runs[K]
        for eps in epsilons:
            tol_e = 1e-12 if eps == 0.0 else eps
            Dm = ref["D_max_seq"] if "D_max_seq" in ref else \
                np.array([r["D_max"] for r in ref["rows"]])
            hit = np.nonzero(Dm <= tol_e)[0]
            t_stop = int(hit[0]) if len(hit) else -1
            if t_stop < 0:
                continue
            for B in sorted(set(B_list) | {K}):
                u = U_by_B[K][B]
                u_init, u_best, u_stop = u[0], u.min(), u[t_stop]
                denom = u_init - u_best
                loss = (u_stop - u_best) / denom if denom > 1e-15 else float("nan")
                stop_rows.append({
                    "K": K,
                    "epsilon": eps,
                    "B": B,
                    "scale_label": "fine (B=K)" if B == K else f"coarse B={B}",
                    "t_stop": t_stop,
                    "U_stop": f"{u_stop:.12e}",
                    "U_best": f"{u_best:.12e}",
                    "t_best": int(u.argmin()),
                    "absolute_gap": f"{u_stop - u_best:.12e}",
                    "residual_regret": f"{(u_stop - u_best) / u_best:.6f}",
                    "improvement_loss": f"{loss:.6f}" if loss == loss else "",
                    "improvement_retention": f"{1.0 - loss:.6f}" if loss == loss else "",
                })
    stop_csv = out_dir / "stop_quality_by_scale.csv"
    write_csv(stop_csv, list(stop_rows[0]), stop_rows)

    # ---------------- 4. regression anchors ----------------
    m1b_csv = REPO_ROOT / "experiments/m1b_tolerance/results/epsilon_k_summary.csv"
    m1b_rows = list(csv.DictReader(m1b_csv.open(encoding="utf-8")))
    m1b_stop = {(int(r["K"]), float(r["epsilon"])): r["stopping_round"]
                for r in m1b_rows}
    ok_stop = True
    for K in K_list:
        Dm = np.array([r["D_max"] for r in runs[K]["rows"]])
        for eps in epsilons:
            tol_e = 1e-12 if eps == 0.0 else eps
            hit = np.nonzero(Dm <= tol_e)[0]
            t_here = int(hit[0]) if len(hit) else "none"
            t_committed = m1b_stop[(K, eps)]
            if str(t_here) != t_committed:
                ok_stop = False
    check_stop = ok_stop
    print(f"\n[anchor] stopping times match committed M1B results: {check_stop}")
    m11_csv = REPO_ROOT / "experiments/m1a_1_long_horizon/results/long_horizon_trajectory_K100.csv"
    m11 = list(csv.DictReader(m11_csv.open(encoding="utf-8")))
    ok_traj = all(
        abs(runs[100]["rows"][t]["U_L2"] - float(m11[t]["U_L2"]))
        <= 1e-12 * abs(float(m11[t]["U_L2"]))
        for t in range(101)
    )
    print(f"[anchor] K=100 trajectory matches committed M1A.1 (t=0..100): {ok_traj}")

    # ---------------- 5. K-stability of fixed-B optima ----------------
    scale_rows: list[dict] = []
    for B in B_list:
        ts = [next(r["t_best"] for r in opt_rows if r["K"] == K and r["B"] == B)
              for K in K_list]
        scale_rows.append({
            "scale": f"coarse B={B}",
            "t_best_K50": ts[0], "t_best_K100": ts[1], "t_best_K200": ts[2],
            "t_best_K400": ts[3],
            "max_min_ratio": f"{max(ts) / min(ts):.3f}",
            "note": "",
        })
    for K in K_list:
        fine_ts = next(r["t_best"] for r in opt_rows if r["K"] == K and r["B"] == K)
        scale_rows.append({
            "scale": f"fine (B=K={K})",
            "t_best_K50": "", "t_best_K100": "", "t_best_K200": "", "t_best_K400": "",
            "max_min_ratio": "",
            "note": f"fine-grid t_best at K={K}: {fine_ts}",
        })
    scale_csv = out_dir / "fixed_B_t_best_vs_K.csv"
    write_csv(scale_csv, list(scale_rows[0]), scale_rows)

    # ---------------- 6. multiscale diagnostic summary (optional, non-primary) ----
    ms_rows: list[dict] = []
    for K in K_list:
        for eps in epsilons:
            tol_e = 1e-12 if eps == 0.0 else eps
            Dm = np.array([r["D_max"] for r in runs[K]["rows"]])
            hit = np.nonzero(Dm <= tol_e)[0]
            t_stop = int(hit[0]) if len(hit) else -1
            if t_stop < 0:
                continue
            rets = [float(r["improvement_retention"]) for r in stop_rows
                    if r["K"] == K and r["epsilon"] == eps and r["B"] in B_list]
            ms_rows.append({
                "K": K, "epsilon": eps,
                "mean_retention_over_B5_50": f"{np.mean(rets):.6f}",
                "min_retention_over_B5_50": f"{np.min(rets):.6f}",
                "max_retention_over_B5_50": f"{np.max(rets):.6f}",
            })
    ms_csv = out_dir / "multiscale_retention_summary.csv"
    write_csv(ms_csv, list(ms_rows[0]), ms_rows)

    # ---------------- 7. figures ----------------
    fig1 = out_dir / "fig_t_best_vs_K.png"
    plot_t_best_vs_K(opt_rows, K_list, B_list, path=fig1)
    fig2 = out_dir / "fig_retention_vs_epsilon.png"
    plot_retention_vs_epsilon(stop_rows, K_list, B_list, epsilons, path=fig2)
    fig3 = out_dir / "fig_regret_vs_loss.png"
    plot_regret_vs_loss(stop_rows, B_list, path=fig3)
    fig4 = out_dir / "fig_coarse_profiles.png"
    st = runs[100]["states_arr"]
    P0 = st[0].reshape(10, 10).sum(axis=1)
    t_e005 = next(int(r["t_stop"]) for r in stop_rows
                  if r["K"] == 100 and float(r["epsilon"]) == 0.005)
    t_best_fine = next(int(r["t_best"]) for r in opt_rows
                       if r["K"] == 100 and r["B"] == 100)
    plot_coarse_profiles(
        {
            "initial t=0": P0,
            "eps=0.005 stop state": st[t_e005].reshape(10, 10).sum(axis=1),
            "fine U-best state": st[t_best_fine].reshape(10, 10).sum(axis=1),
        },
        B=10, path=fig4,
    )

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
            "anchors_pass": {"stopping_times_match_M1B": check_stop,
                             "K100_matches_M1A1": ok_traj},
            "fixed_B_t_best_stable": {
                str(r["scale"]): r.get("max_min_ratio", "")
                for r in scale_rows if str(r["scale"]).startswith("coarse")
            },
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m1b_1_metric_scale/run_experiment.py",
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
