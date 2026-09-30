"""Runner for M1C.5: boundary-gated rebound theory verification.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m1c5_boundary_gate/run_analysis.py

Verification of the M1C.5 analytic boundary gate on the EXISTING canonical
Variant A trajectories (K in {50,100,200,400}; no new model experiment):

  * Region-1 reduction identity: M*Delta - gamma == (1-alpha)D - D_max
    + alpha*h  (checked per round, max residual reported);
  * state-reduced necessary gate: every rebound round must satisfy
    a_t > a_crit(D_max,t) with a_crit(D) = [(1-D)+sqrt(D(D+2))]/2
    (0 false negatives expected — theorem);
  * gate-open rounds (a > a_crit) that still contract — necessary is not
    sufficient; their profile factor P(x_r) = D(x_r)/D_max is compared
    with rebound rounds;
  * cross-K a_crit quantiles within the churn window (M1C.4 definition).
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

from sand_m0.adaptive import cumulative_excess, run_m1a_deterministic  # noqa: E402
from sand_m0.kernel_theory import a_crit_gate, state_reduced_gate  # noqa: E402
from sand_m0.model import DISTRIBUTIONS  # noqa: E402
from sand_m0.plotting import (  # noqa: E402
    plot_gate_scatter,
    plot_gate_margin,
    plot_profile_factor,
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
    parser.add_argument("--config", type=Path,
                        default=Path(__file__).resolve().parent / "config.json")
    args = parser.parse_args(argv)
    config_path = args.config.resolve()
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    out_dir = config_path.parent / cfg["output_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    K_list = list(cfg["K_list"])
    alpha, tol = cfg["alpha"], cfg["stop_tolerance_exact"]
    q_near = DISTRIBUTIONS["near"]
    fp_tol = cfg["fp_tolerance"]
    started = time.time()

    print(f"M1C.5 boundary-gate verification — config: {config_path}")
    expected_stops = {50: 324, 100: 699, 200: 1445, 400: 2876}

    per_rows: list[dict] = []
    summaries: list[dict] = []
    worst_resid = 0.0
    worst_hgap = 0.0
    fn_total = 0
    for K in K_list:
        run = run_m1a_deterministic(q_near, K=K, T_max=cfg["T_max"], alpha=alpha,
                                    stop_tolerance=tol)
        assert run["stopped_at"] == expected_stops[K], "baseline drift"
        stop = run["stopped_at"]
        states = np.asarray(run["states"], dtype=float)
        rows = run["rows"]
        # churn window onset (M1C.4/M1B.2 convention)
        Dm_all = np.array([r["D_max"] for r in rows[:stop]])
        onset = int(np.nonzero(Dm_all < 0.05)[0][0])

        recs = []
        for t in range(stop):
            mass = states[t]
            D = cumulative_excess(mass, K)
            d_max = float(D.max())
            j_star = int(np.argmax(D)) + 1
            a = j_star / K
            gate = state_reduced_gate(mass, j_star, alpha)
            # residual + realized driving prefix (argmax of M*Delta - gamma)
            gamma = d_max - _tm(mass, j_star, alpha)
            g_bins = q_near.bin_probs(K)
            delta = np.cumsum(g_bins) - np.arange(1, K + 1) / K
            moved = alpha * gate["F_a"]
            resid = moved * delta - gamma
            worst_resid = max(worst_resid,
                              float(np.abs(resid[:j_star] - _reduced(mass, j_star, alpha, g_bins)[:j_star]).max()))
            h_gap = abs(gate["h_grid_max"] - gate["h_closed_form"]) \
                if gate["x_star"] is not None else 0.0
            worst_hgap = max(worst_hgap, h_gap)
            j_r = int(np.argmax(resid)) + 1
            rebound = bool(resid.max() > 0.0)
            B = a - gate["a_crit"]
            recs.append({
                "K": K, "t": t, "D_max": f"{d_max:.10e}",
                "a_t": f"{a:.5f}", "F_a": f"{gate['F_a']:.8f}",
                "a_crit": f"{gate['a_crit']:.6f}", "B_t": f"{B:.6f}",
                "gate_open": B > 0.0, "rebound": rebound,
                "x_star": (f"{gate['x_star']:.5f}" if gate["x_star"] is not None else ""),
                "h_grid_max": f"{gate['h_grid_max']:.8f}",
                "h_closed_form": (f"{gate['h_closed_form']:.8f}"
                                  if gate["x_star"] is not None else ""),
                "x_r": f"{j_r / K:.5f}", "j_r_le_jstar": j_r <= j_star,
                "P_profile": (f"{float(D[j_r - 1]) / d_max:.6f}" if d_max > 0 else ""),
                "residual_max": f"{float(resid.max()):.8e}",
                "in_churn_window": t >= onset,
            })
            if rebound and B <= fp_tol:
                fn_total += 1  # theorem says this must stay 0
        per_rows.extend(recs)

        # summary over the churn window
        win = recs[onset:]
        ac = np.array([float(r["a_crit"]) for r in win])
        B = np.array([float(r["B_t"]) for r in win])
        reb = np.array([r["rebound"] for r in win], dtype=bool)
        gate_open = np.array([r["gate_open"] for r in win], dtype=bool)
        P_reb = np.array([float(r["P_profile"]) for r in win
                          if r["rebound"] and r["P_profile"] != ""])
        P_open_con = np.array([float(r["P_profile"]) for r in win
                               if r["gate_open"] and not r["rebound"]
                               and r["P_profile"] != ""])
        summaries.append({
            "K": K, "t_stop": stop, "t_onset_0.05": onset,
            "window_rounds": len(win),
            "acrit_p5": f"{float(np.percentile(ac, 5)):.5f}",
            "acrit_p25": f"{float(np.percentile(ac, 25)):.5f}",
            "acrit_p50": f"{float(np.percentile(ac, 50)):.5f}",
            "acrit_p75": f"{float(np.percentile(ac, 75)):.5f}",
            "acrit_p95": f"{float(np.percentile(ac, 95)):.5f}",
            "rebound_count": int(reb.sum()),
            "rebound_fraction": f"{float(reb.mean()):.5f}",
            "false_negatives": int((reb & ~gate_open).sum()),
            "gate_open_rounds": int(gate_open.sum()),
            "gate_open_contraction": int((gate_open & ~reb).sum()),
            "gate_open_rebound": int((gate_open & reb).sum()),
            "gate_closed_rounds": int((~gate_open).sum()),
            "P_rebound_median": f"{float(np.median(P_reb)):.5f}" if len(P_reb) else "",
            "P_gateopen_contraction_median":
                f"{float(np.median(P_open_con)):.5f}" if len(P_open_con) else "",
        })
        print(f"  [K={K}] stop={stop} window={len(win)} rebounds={int(reb.sum())} "
              f"a_crit p50={summaries[-1]['acrit_p50']} FN={summaries[-1]['false_negatives']} "
              f"gate_open_contraction={summaries[-1]['gate_open_contraction']}")

    # ---------------- CSV outputs ----------------
    write_csv(out_dir / "per_round_gate.csv", list(per_rows[0]), per_rows)
    write_csv(out_dir / "gate_summary.csv", list(summaries[0]), summaries)
    print("\n=== boundary gate summary (churn window) ===")
    for s in summaries:
        print(f"  K={s['K']:>3}: a_crit p5/25/50/75/95 = {s['acrit_p5']}/{s['acrit_p25']}/"
              f"{s['acrit_p50']}/{s['acrit_p75']}/{s['acrit_p95']}  "
              f"reb={float(s['rebound_fraction']):.4f}  FN={s['false_negatives']}  "
              f"gate_open={s['gate_open_rounds']} (reb {s['gate_open_rebound']} / "
              f"con {s['gate_open_contraction']})  "
              f"P_reb={s['P_rebound_median']} P_open_con={s['P_gateopen_contraction_median']}")
    print(f"\nworst Region-1 reduction residual: {worst_resid:.2e}")
    print(f"worst |h_grid_max - h_closed_form|: {worst_hgap:.2e} (O(1/K^2) rounding)")
    print(f"total false negatives of the necessary gate: {fn_total} (theorem: 0)")

    # ---------------- figures ----------------
    plot_gate_scatter(per_rows, K_panels=[100, 400], path=out_dir / "fig_gate_scatter.png")
    plot_gate_margin(per_rows, path=out_dir / "fig_gate_margin_B.png")
    plot_profile_factor(per_rows, path=out_dir / "fig_profile_factor.png")

    # ---------------- metadata ----------------
    output_files = sorted(p.name for p in out_dir.iterdir()
                          if p.is_file() and p.name != "metadata.json")
    metadata = {
        "model": cfg["model"], "model_version": cfg["model_version"],
        "experiment": cfg["experiment"],
        "run_date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": cfg,
        "headline": {
            "worst_region1_reduction_residual": worst_resid,
            "worst_h_grid_vs_closed_form": worst_hgap,
            "necessary_gate_false_negatives_total": fn_total,
            "acrit_p50_by_K": {str(s["K"]): s["acrit_p50"] for s in summaries},
        },
        "environment": {"python": sys.version, "platform": platform.platform(),
                        "numpy": np.__version__},
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m1c5_boundary_gate/run_analysis.py",
        "outputs": {name: sha256_of(out_dir / name) for name in output_files},
    }
    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\noutputs written to: {out_dir} ({len(output_files) + 1} files)")
    print(f"runtime: {metadata['runtime_seconds']} s")
    return 0


def _tm(mass, j_star, alpha):
    from sand_m0.kernel_theory import target_matched_excess
    return target_matched_excess(mass, j_star, alpha)


def _reduced(mass, j_star, alpha, g_bins):
    """Region-1 reduced residual form (1-alpha)D - D_max + alpha*h."""
    K = mass.shape[0]
    D = cumulative_excess(mass, K)
    d_max = float(D.max())
    F_a = float(mass[:j_star].sum())
    x = np.arange(1, K + 1) / K
    h = x * (F_a * (2.0 - x) - 1.0)
    return (1.0 - alpha) * D - d_max + alpha * h


if __name__ == "__main__":
    raise SystemExit(main())
