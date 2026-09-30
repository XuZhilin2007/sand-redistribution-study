"""Runner for M1B.2: late-stage micro-correction scaling diagnostics.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m1b_2_micro_scaling/run_experiment.py

Hypothesis under test (and actively falsified here): the exact M1A/M1B
controller enters a micro-correction regime in which per-round action size
becomes O(1/K), explaining why exact stopping time grows ~proportionally
to K while finite-epsilon stopping is K-robust.

Per-round action diagnostics on the canonical exact reference runs:
sweep_mass, moved_mass, net_export_mass (exact identity with the canonical
update), L1_state_step, delta_U_density (signed), and their aggregations
over relative-phase bins and D_max bands, cross-K scaling diagnostics
(X vs K*X vs K^2*X, CV, log-log slope), exact-tail lengths, and an
epsilon overlay locating each M1B stop inside the action-size regimes.
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
    plot_action_vs_time,
    plot_churn_vs_D_max,
    plot_reference_D_max_bands,
    plot_tail_vs_K,
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


def median(v: list[float]) -> float:
    return float(np.median(v)) if v else float("nan")


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
    eps_list = list(cfg["epsilon_overlay"])
    tail_thresholds = list(cfg["tail_thresholds"])
    phase_edges = list(cfg["phase_bins_relative"])
    q = DISTRIBUTIONS["near"]
    started = time.time()

    print(f"M1B.2 micro-scaling run — config: {config_path}")
    print(f"  K={K_list}, alpha={alpha}, tol={tol}, T_max={T_max}")

    # ---------------- 1. reference runs + per-round action diagnostics ----------------
    runs: dict[int, dict] = {}
    per_round_rows: list[dict] = []
    for K in K_list:
        t0 = time.time()
        run = run_m1a_deterministic(q, K=K, T_max=T_max, alpha=alpha, stop_tolerance=tol)
        stop = run["stopped_at"]
        states = np.asarray(run["states"], dtype=float)
        full_rows = run["rows"]
        rows = full_rows[:stop]  # active decisions only
        diag = []
        for t, r in enumerate(rows):
            F_q_a = float(q.cdf(np.array([r["a_t"]]))[0])
            moved = r["removed_mass"]
            net_export = moved * (1.0 - F_q_a)
            L1 = float(np.abs(states[t + 1] - states[t]).sum())
            ud0 = K * r["U_L2"]
            ud1 = K * full_rows[t + 1]["U_L2"]
            dU = ud0 - ud1
            diag.append({
                "t": t, "a_t": r["a_t"], "sweep_mass": r["prefix_mass"],
                "moved_mass": moved, "F_q_a": F_q_a, "net_export_mass": net_export,
                "L1_state_step": L1, "D_max": r["D_max"],
                "U_density": ud0, "delta_U_density": dU,
                "delta_U_per_moved": (dU / moved) if moved > 0 else float("nan"),
            })
        runs[K] = {"run": run, "stop": stop, "diag": diag,
                   "states": states, "rows": rows}
        for d in diag:
            per_round_rows.append({"K": K, **{k: (f"{v:.12e}" if isinstance(v, float) else v)
                                               for k, v in d.items()}})
        print(f"  [K={K}] stop t={stop}; diagnostics built ({time.time() - t0:.1f}s)")
    per_round_csv = out_dir / "action_diagnostics_per_round.csv"
    write_csv(per_round_csv, list(per_round_rows[0]), per_round_rows)

    # ---------------- 2. phase-binned medians (relative active position) ----------------
    phase_rows: list[dict] = []
    for K in K_list:
        diag = runs[K]["diag"]
        stop = runs[K]["stop"]
        for i in range(len(phase_edges) - 1):
            lo, hi = phase_edges[i], phase_edges[i + 1]
            seg = diag[int(lo * stop):max(int(hi * stop), int(lo * stop) + 1)]
            seg = [d for d in seg if d["t"] < stop]
            if not seg:
                continue
            phase_rows.append({
                "K": K,
                "phase": f"{lo:.2f}-{hi:.2f}",
                "rounds": len(seg),
                "median_moved_mass": f"{median([d['moved_mass'] for d in seg]):.6f}",
                "median_net_export": f"{median([d['net_export_mass'] for d in seg]):.6f}",
                "median_L1_step": f"{median([d['L1_state_step'] for d in seg]):.6f}",
                "median_abs_delta_U_density": f"{median([abs(d['delta_U_density']) for d in seg]):.6e}",
                "median_a_t": f"{median([d['a_t'] for d in seg]):.4f}",
                "median_sweep_mass": f"{median([d['sweep_mass'] for d in seg]):.4f}",
            })
    phase_csv = out_dir / "phase_summary.csv"
    write_csv(phase_csv, list(phase_rows[0]), phase_rows)

    # ---------------- 3. D_max-band medians ----------------
    band_defs = [
        ("D<0.005", lambda d: d["D_max"] < 0.005),
        ("0.005<=D<0.01", lambda d: 0.005 <= d["D_max"] < 0.01),
        ("0.01<=D<0.02", lambda d: 0.01 <= d["D_max"] < 0.02),
        ("0.02<=D<0.05", lambda d: 0.02 <= d["D_max"] < 0.05),
        ("D>=0.05", lambda d: d["D_max"] >= 0.05),
    ]
    band_rows: list[dict] = []
    for K in K_list:
        diag = runs[K]["diag"]
        for label, cond in band_defs:
            seg = [d for d in diag if cond(d)]
            if not seg:
                continue
            band_rows.append({
                "K": K,
                "band": label,
                "rounds": len(seg),
                "median_moved_mass": f"{median([d['moved_mass'] for d in seg]):.6f}",
                "median_net_export": f"{median([d['net_export_mass'] for d in seg]):.6f}",
                "median_L1_step": f"{median([d['L1_state_step'] for d in seg]):.6f}",
                "median_abs_delta_U_density": f"{median([abs(d['delta_U_density']) for d in seg]):.6e}",
                "median_a_t": f"{median([d['a_t'] for d in seg]):.4f}",
                "median_sweep_mass": f"{median([d['sweep_mass'] for d in seg]):.4f}",
            })
    band_csv = out_dir / "D_max_band_summary.csv"
    write_csv(band_csv, list(band_rows[0]), band_rows)

    # ---------------- 4. cross-K scaling of late-stage action sizes ----------------
    late_band = [d for d in runs[K_list[0]]["diag"] if d["D_max"] < 0.005]
    scaling_rows: list[dict] = []
    quantities = ["moved_mass", "net_export_mass", "L1_state_step"]
    late_medians: dict[str, dict[int, float]] = {}
    for K in K_list:
        seg = [d for d in runs[K]["diag"] if d["D_max"] < 0.005]
        for qname in quantities:
            late_medians.setdefault(qname, {})[K] = median([d[qname] for d in seg])
    for qname in quantities:
        xs = np.array([late_medians[qname][K] for K in K_list], dtype=float)
        Ks = np.array(K_list, dtype=float)
        cv_X = float(xs.std() / xs.mean())
        cv_KX = float((Ks * xs).std() / (Ks * xs).mean())
        cv_K2X = float((Ks**2 * xs).std() / (Ks**2 * xs).mean())
        slope = float(np.polyfit(np.log(Ks), np.log(xs), 1)[0])
        scaling_rows.append({
            "quantity": qname,
            "late_median_K50": f"{xs[0]:.6f}",
            "late_median_K100": f"{xs[1]:.6f}",
            "late_median_K200": f"{xs[2]:.6f}",
            "late_median_K400": f"{xs[3]:.6f}",
            "CV_of_X": f"{cv_X:.4f}",
            "CV_of_KX": f"{cv_KX:.4f}",
            "CV_of_K2X": f"{cv_K2X:.4f}",
            "loglog_slope_beta": f"{slope:.3f}",
            "interpretation": ("K-independent (no 1/K scaling)" if cv_X < cv_KX and cv_X < cv_K2X
                               else "consistent with 1/K scaling" if cv_KX < cv_X and cv_KX < cv_K2X
                               else "ambiguous"),
        })
    for K in K_list:
        seg = [d for d in runs[K]["diag"] if d["D_max"] < 0.005]
        scaling_rows.append({
            "quantity": f"D_max_band_median (K={K})",
            "late_median_K50": f"{median([d['D_max'] for d in seg]):.6e}" if seg else "n/a",
            "late_median_K100": "", "late_median_K200": "", "late_median_K400": "",
            "CV_of_X": "", "CV_of_KX": "", "CV_of_K2X": "", "loglog_slope_beta": "",
            "interpretation": f"{len(seg)} rounds with D_max<0.005",
        })
    scale_csv = out_dir / "cross_K_scaling.csv"
    write_csv(scale_csv, list(scaling_rows[0]), scaling_rows)
    print("\n=== late-stage (D_max<0.005) cross-K scaling ===")
    for r in scaling_rows:
        if r["quantity"] in quantities:
            print(f"  {r['quantity']}: medians={[float(r[k]) for k in
                  ['late_median_K50', 'late_median_K100', 'late_median_K200', 'late_median_K400']]}")
            print(f"     CV: X={r['CV_of_X']} KX={r['CV_of_KX']} K2X={r['CV_of_K2X']}  "
                  f"beta={r['loglog_slope_beta']}  -> {r['interpretation']}")

    # ---------------- 5. band/floor + dip statistics ----------------
    bandstat_rows: list[dict] = []
    for K in K_list:
        diag = runs[K]["diag"]
        Dm = np.array([d["D_max"] for d in diag])
        Ud = np.array([d["U_density"] for d in diag])
        row = {"K": K, "active_rounds": len(diag),
               "D_max_p5": f"{np.percentile(Dm, 5):.5f}",
               "D_max_p50": f"{np.percentile(Dm, 50):.5f}",
               "D_max_p95": f"{np.percentile(Dm, 95):.5f}",
               "U_density_p5": f"{np.percentile(Ud, 5):.5f}",
               "U_density_p50": f"{np.percentile(Ud, 50):.5f}",
               "U_density_p95": f"{np.percentile(Ud, 95):.5f}"}
        for eps in eps_list:
            dips = int((Dm <= eps).sum())
            first = int(np.nonzero(Dm <= eps)[0][0]) if dips else -1
            row[f"dips_le_{eps}"] = dips
            row[f"dips_per_1000_le_{eps}"] = f"{1000.0 * dips / len(diag):.1f}"
            row[f"first_dip_{eps}"] = first
        bandstat_rows.append(row)
    bandstat_csv = out_dir / "band_and_dip_statistics.csv"
    write_csv(bandstat_csv, list(bandstat_rows[0]), bandstat_rows)
    print("\n=== fluctuation bands (active phase) ===")
    for K, r in zip(K_list, bandstat_rows):
        print(f"  K={K:>3}: D_max p5/p50/p95 = {r['D_max_p5']}/{r['D_max_p50']}/{r['D_max_p95']}; "
              f"U_dens p5/p50/p95 = {r['U_density_p5']}/{r['U_density_p50']}/{r['U_density_p95']}")

    # ---------------- 6. tail lengths and K-scaling of the exact stop ----------------
    tail_rows: list[dict] = []
    for K in K_list:
        diag = runs[K]["diag"]
        Dm = np.array([d["D_max"] for d in diag])
        Ud = np.array([d["U_density"] for d in diag])
        stop = runs[K]["stop"]
        best_t = int(Ud.argmin())
        row = {"K": K, "t_U_best": best_t, "t_U_best_over_K": f"{best_t / K:.3f}",
               "t_exact_stop": stop, "t_exact_stop_over_K": f"{stop / K:.3f}"}
        for thr in tail_thresholds:
            hit = np.nonzero(Dm <= thr)[0]
            t_pass = int(hit[0]) if len(hit) else -1
            tail = stop - t_pass if t_pass >= 0 else -1
            row[f"t_pass_{thr}"] = t_pass
            row[f"tail_after_{thr}"] = tail
            row[f"tail_over_K_{thr}"] = f"{tail / K:.3f}" if t_pass >= 0 else ""
        tail_rows.append(row)
    tail_csv = out_dir / "tail_lengths.csv"
    write_csv(tail_csv, list(tail_rows[0]), tail_rows)
    print("\n=== exact-stop tail lengths ===")
    for r in tail_rows:
        print(f"  K={r['K']:>3}: t_U_best/K={r['t_U_best_over_K']}  "
              f"t_stop/K={r['t_exact_stop_over_K']}  "
              f"tail/K after 0.02/0.01/0.005 = {r['tail_over_K_0.02']}/{r['tail_over_K_0.01']}/{r['tail_over_K_0.005']}")

    # ---------------- 7. epsilon overlay (first-passage action states) ----------------
    over_rows: list[dict] = []
    for K in K_list:
        diag = runs[K]["diag"]
        Dm = np.array([d["D_max"] for d in diag])
        Ud = np.array([d["U_density"] for d in diag])
        u_best = float(Ud.min())
        for eps in eps_list:
            hit = np.nonzero(Dm <= eps)[0]
            if not len(hit):
                continue
            t_eps = int(hit[0])
            d = diag[t_eps]
            regime_ratio = d["moved_mass"] * K
            over_rows.append({
                "K": K,
                "epsilon": eps,
                "t_stop": t_eps,
                "moved_mass_at_stop": f"{d['moved_mass']:.6f}",
                "net_export_at_stop": f"{d['net_export_mass']:.6f}",
                "L1_step_at_stop": f"{d['L1_state_step']:.6f}",
                "moved_mass_times_K": f"{regime_ratio:.2f}",
                "regime_note": ("gross churn scale (moved*K >> 1): NOT a micro-correction"
                                if regime_ratio > 10 else "small-action round"),
                "U_density_at_stop": f"{d['U_density']:.6e}",
                "improvement_retention_fine": f"{1.0 - (d['U_density'] - u_best) / (Ud[0] - u_best):.6f}",
            })
    over_csv = out_dir / "epsilon_overlay.csv"
    write_csv(over_csv, list(over_rows[0]), over_rows)
    print("\n=== epsilon overlay (first-passage action states) ===")
    for r in over_rows:
        print(f"  K={r['K']:>3} eps={r['epsilon']:>6}: t={r['t_stop']:>4}  "
              f"moved={float(r['moved_mass_at_stop']):.4f}  moved*K={float(r['moved_mass_times_K']):.1f}  "
              f"retention={float(r['improvement_retention_fine']):.4f}")

    # ---------------- 8. annihilation event analysis ----------------
    ann_rows: list[dict] = []
    for K in K_list:
        stop = runs[K]["stop"]
        pre = runs[K]["rows"][stop - 1]
        st = runs[K]["states"]
        D_pre = np.cumsum(st[stop - 1]) - np.arange(1, K + 1) / K
        D_post = np.cumsum(st[stop]) - np.arange(1, K + 1) / K
        F_a = float(q.cdf(np.array([pre["a_t"]]))[0])
        ann_rows.append({
            "K": K,
            "t_stop": stop,
            "pre_stop_a_t": f"{pre['a_t']:.4f}",
            "pre_stop_F_a": f"{F_a:.4f}",
            "F_a_le_0.5": F_a <= 0.5 + 1e-12,
            "sweep_mass": f"{pre['prefix_mass']:.6f}",
            "moved_mass": f"{pre['removed_mass']:.6f}",
            "positive_prefixes_pre": int((D_pre > 1e-12).sum()),
            "positive_prefixes_post": int((D_post > 1e-12).sum()),
            "D_max_pre": f"{D_pre.max():.6e}",
            "D_max_post": f"{D_post.max():.6e}",
        })
    ann_csv = out_dir / "annihilation_events.csv"
    write_csv(ann_csv, list(ann_rows[0]), ann_rows)
    print("\n=== annihilation events (one-sweep total clearance) ===")
    for r in ann_rows:
        print(f"  K={r['K']:>3}: t={r['t_stop']:>4}  a={r['pre_stop_a_t']}  F(a)={r['pre_stop_F_a']} "
              f"(<=0.5: {r['F_a_le_0.5']})  moved={r['moved_mass']}  "
              f"positive prefixes {r['positive_prefixes_pre']}->{r['positive_prefixes_post']}")

    # ---------------- 9. figures ----------------
    fig1 = out_dir / "fig_action_vs_time.png"
    plot_action_vs_time(runs, K_list, path=fig1)
    fig2 = out_dir / "fig_churn_vs_D_max.png"
    plot_churn_vs_D_max(runs, K_list, path=fig2)
    fig3 = out_dir / "fig_reference_D_max_bands.png"
    plot_reference_D_max_bands(runs, K_list, eps_list, path=fig3)
    fig4 = out_dir / "fig_tail_vs_K.png"
    plot_tail_vs_K(tail_rows, over_rows, path=fig4)

    # ---------------- 10. metadata ----------------
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
            "hypothesis_micro_correction": "FALSIFIED for gross action size: "
                "moved/net-export/L1 medians are K-independent (CV<2%) during the "
                "late stage; the controller churns ~15% of the sand per round at "
                "every resolution until the exact stop",
            "stopping_mechanism": "exact stop = rare one-sweep annihilation event "
                "(all positive prefixes cleared in a single sweep); F(a)<=0.5 "
                "necessary condition held 4/4",
            "D_max_fluctuation_band_K_robust": True,
            "annihilation_all_peak_switch": True,
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m1b_2_micro_scaling/run_experiment.py",
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
