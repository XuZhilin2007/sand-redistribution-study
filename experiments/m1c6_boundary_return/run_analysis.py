"""Runner for M1C.6: boundary return & peak-competition dynamics.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m1c6_boundary_return/run_analysis.py

Audit of the existing canonical Variant A deterministic trajectories
(K in {50,100,200,400}; no new model). Closes the churn feedback loop:

  rebound (gate open) -> peak/boundary switch -> contractions
  -> boundary returns across the analytic gate -> next rebound becomes possible.

Analyses:
  * rebound episodes (rebound -> next rebound; terminal episode -> stop);
  * gate-exit/return: B_t = a_t - a_crit(D_max,t) aligned at each rebound,
    tau_return / tau_rebound distributions, empirical episode frequencies;
  * peak competition: meaningful local maxima (plateau-safe, sep = max(2, K//50)),
    the exact pairwise peak-height update identity (M1C.2 corollary),
    active-vs-competing peak motion in contraction rounds;
  * rebound-induced switch mechanism: P1 (existing secondary overtakes) vs
    P_multi (third peak) vs P2 (reinjection-rescaled new peak) vs no switch;
  * boundary map a_t -> a_{t+1} by round type;
  * terminal episode vs ordinary episodes.
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
from sand_m0.kernel_theory import a_crit_gate, target_matched_excess  # noqa: E402
from sand_m0.model import DISTRIBUTIONS  # noqa: E402
from sand_m0.plotting import (  # noqa: E402
    plot_gate_episode_alignment,
    plot_boundary_map,
    plot_peak_gap_alignment,
    plot_terminal_vs_ordinary,
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


def local_peaks(mass: np.ndarray):
    """Local maxima of the cumulative-excess profile (plateau -> first index)."""
    D = cumulative_excess(mass, mass.shape[0])
    pk: list[int] = []
    for j in range(mass.shape[0]):
        left = D[j - 1] if j > 0 else -np.inf
        right = D[j + 1] if j < mass.shape[0] - 1 else -np.inf
        if D[j] >= right and (j == 0 or D[j] > left):
            if pk and j == pk[-1] + 1 and D[j] == D[pk[-1]]:
                continue
            pk.append(j)
    return D, pk


def analyze_run(K: int, run: dict, alpha: float, g_bins: np.ndarray,
                cfg: dict) -> dict:
    stop = run["stopped_at"]
    states = np.asarray(run["states"], dtype=float)
    sep = max(2, K // 50)
    sw_thr = cfg["large_switch_threshold"]
    delta = np.cumsum(g_bins) - np.arange(1, K + 1) / K

    # per-round arrays
    a = np.empty(stop + 1)
    Dm = np.empty(stop + 1)
    peaks = []
    for t in range(stop + 1):
        D, pk = local_peaks(states[t])
        Dm[t] = float(D.max())
        j1 = pk[int(np.argmax([D[j] for j in pk]))] if pk else int(np.argmax(D))
        a[t] = (j1 + 1) / K
        far = [j for j in pk if abs(j - j1) > sep]
        j2 = far[int(np.argmax([D[j] for j in far]))] if far else None
        peaks.append({"j1": j1, "D1": float(D[j1]),
                      "j2": j2, "D2": (float(D[j2]) if j2 is not None else None),
                      "n_pk": len(pk)})
    # clamp tiny negative D_max (stopped states) to 0: the gate is
    # inert there and a_crit(0) = 1/2
    B = a - np.array([a_crit_gate(max(d, 0.0)) for d in Dm])
    reb = np.array([bool(Dm[t + 1] > Dm[t] + cfg["fp_tolerance"])
                    for t in range(stop)])
    t_r_list = [t for t in range(stop) if reb[t]]

    # episode segmentation + gate return
    episodes = []
    for ep, t_r in enumerate(t_r_list):
        nxt = t_r_list[ep + 1] if ep + 1 < len(t_r_list) else stop
        terminal = ep + 1 >= len(t_r_list)
        tau_return = -1
        for s in range(1, nxt - t_r + 1):
            if B[t_r + s] > 0:
                tau_return = s
                break
        # switch mechanism (pre-rebound peak structure of state t_r)
        pk_info = peaks[t_r]
        j1, j2 = pk_info["j1"], pk_info["j2"]
        a_pre, a_post = a[t_r], a[t_r + 1]
        j_post = int(np.argmax(cumulative_excess(states[t_r + 1], K)))
        D_pre, pk_all = local_peaks(states[t_r])
        mech = "no_switch"
        if abs(a_post - a_pre) > sw_thr:
            if j2 is not None and abs(j_post - j2) <= 2:
                mech = "P1"
            elif any(abs(j_post - j) <= 2 for j in pk_all
                     if abs(j - j1) > sep and j != j2):
                mech = "P_multi"
            else:
                mech = "P2_new"
        episodes.append({
            "K": K, "episode": ep, "t_r": t_r, "terminal": terminal,
            "a_r": f"{a[t_r]:.5f}", "B_r": f"{B[t_r]:.6f}",
            "a_post": f"{a[t_r + 1]:.5f}", "B_post": f"{B[t_r + 1]:.6f}",
            "switch": abs(a_post - a_pre) > sw_thr,
            "mechanism": mech,
            "x2_pre": (f"{(j2 + 1) / K:.5f}" if j2 is not None else ""),
            "H_pre": (f"{pk_info['D2'] - pk_info['D1']:.6f}"
                      if j2 is not None else ""),
            "tau_return": tau_return,
            "tau_rebound": nxt - t_r,
            "D_max_r": f"{Dm[t_r]:.8f}",
            "in_churn_window": t_r >= 23,
        })

    # pairwise peak update identity + contraction-round peak motion
    worst_pair_resid = 0.0
    d_act_l, d_comp_l = [], []
    for t in range(stop):
        pk = peaks[t]
        j1, j2 = pk["j1"], pk["j2"]
        if j2 is None:
            continue
        D = cumulative_excess(states[t], K)
        Dn = cumulative_excess(states[t + 1], K)
        M = alpha * float(states[t][:j1 + 1].sum())  # 1-based prefix j_star = j1+1
        # kernel j_star is 1-based; j1/j2 from local_peaks are 0-based
        tm_prof = target_matched_excess(states[t], j1 + 1, alpha)
        tm1, tm2 = tm_prof[j1], tm_prof[j2]
        pred = (tm1 - tm2) + M * (delta[j1] - delta[j2])
        worst_pair_resid = max(worst_pair_resid,
                               abs(float(Dn[j1] - Dn[j2]) - pred))
        if not reb[t]:
            d_act_l.append(float(Dn[j1] - D[j1]))
            d_comp_l.append(float(Dn[j2] - D[j2]))

    # alignment of B at rebounds
    align = {}
    for s in range(cfg["alignment_window"][0], cfg["alignment_window"][1] + 1):
        vals = [B[t_r + s] for t_r in t_r_list if 0 <= t_r + s <= stop]
        align[s] = np.array(vals)

    # boundary map
    map_rows = []
    for t in range(stop):
        map_rows.append({
            "K": K, "t": t, "a_t": f"{a[t]:.5f}", "a_next": f"{a[t + 1]:.5f}",
            "da": f"{a[t + 1] - a[t]:+.5f}", "rebound": bool(reb[t]),
            "gate_open_t": B[t] > 0, "B_t": f"{B[t]:.6f}",
            "in_churn_window": t >= 23,
        })

    return {
        "stop": stop, "a": a, "Dm": Dm, "B": B, "reb": reb,
        "t_r_list": t_r_list, "episodes": episodes,
        "peaks": peaks, "align": align, "map_rows": map_rows,
        "worst_pair_resid": worst_pair_resid,
        "d_act": np.array(d_act_l), "d_comp": np.array(d_comp_l),
        "onset": int(np.nonzero(Dm[:stop] < 0.05)[0][0]),
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
    started = time.time()
    expected_stops = {50: 324, 100: 699, 200: 1445, 400: 2876}

    print(f"M1C.6 boundary-return audit — config: {config_path}")
    res = {}
    for K in K_list:
        run = run_m1a_deterministic(q_near, K=K, T_max=cfg["T_max"], alpha=alpha,
                                    stop_tolerance=tol)
        assert run["stopped_at"] == expected_stops[K], "baseline drift"
        res[K] = analyze_run(K, run, alpha, q_near.bin_probs(K), cfg)
        print(f"  [K={K}] stop={res[K]['stop']} episodes={len(res[K]['t_r_list'])} "
              f"pair-identity residual={res[K]['worst_pair_resid']:.2e}")

    # ---- episode CSVs ----
    ep_rows = [e for K in K_list for e in res[K]["episodes"]]
    write_csv(out_dir / "episode_table.csv", list(ep_rows[0]), ep_rows)
    write_csv(out_dir / "boundary_map.csv", list(res[K_list[0]]["map_rows"][0]),
              [r for K in K_list for r in res[K]["map_rows"]])

    # ---- gate-return summary ----
    sum_rows = []
    for K in K_list:
        v = res[K]
        eps = [e for e in v["episodes"] if not e["terminal"]]
        taus = [e["tau_return"] for e in eps]
        taus_pos = [t for t in taus if t > 0]
        mech = {"P1": 0, "P_multi": 0, "P2_new": 0, "no_switch": 0}
        for e in eps:
            mech[e["mechanism"]] += 1
        n = len(eps)
        tr = np.array([e["tau_rebound"] for e in eps])
        sum_rows.append({
            "K": K, "episodes": n,
            "P_post_B_neg": f"{np.mean([float(e['B_post']) < 0 for e in eps]):.4f}",
            "tau_return_median": f"{np.median(taus_pos):.2f}" if taus_pos else "",
            "return_in_1": f"{np.mean([t == 1 for t in taus]):.4f}",
            "return_in_2": f"{np.mean([t == 2 for t in taus]):.4f}",
            "return_in_3p": f"{np.mean([t >= 3 for t in taus]):.4f}",
            "return_never": f"{np.mean([t == -1 for t in taus]):.4f}",
            "tau_rebound_p25": f"{np.percentile(tr, 25):.2f}",
            "tau_rebound_median": f"{np.median(tr):.2f}",
            "tau_rebound_p75": f"{np.percentile(tr, 75):.2f}",
            "tau_rebound_max": int(tr.max()),
            "frac_P1": f"{mech['P1'] / n:.4f}", "frac_P_multi": f"{mech['P_multi'] / n:.4f}",
            "frac_P2_new": f"{mech['P2_new'] / n:.4f}",
            "frac_no_switch": f"{mech['no_switch'] / n:.4f}",
            "pair_identity_residual": f"{v['worst_pair_resid']:.2e}",
            "gap_closing_share_contraction":
                f"{np.mean(v['d_comp'] - v['d_act'] > 0):.4f}",
            "d_active_median": f"{np.median(v['d_act']):+.6f}",
            "d_competing_median": f"{np.median(v['d_comp']):+.6f}",
        })
    write_csv(out_dir / "gate_return_summary.csv", list(sum_rows[0]), sum_rows)
    print("\n=== gate-return summary (episodes) ===")
    for s in sum_rows:
        print(f"  K={s['K']:>3}: eps={s['episodes']:>4}  P(B_post<0)={s['P_post_B_neg']}  "
              f"return 1/2/3+/never = {s['return_in_1']}/{s['return_in_2']}/"
              f"{s['return_in_3p']}/{s['return_never']}  "
              f"tau_reb med={s['tau_rebound_median']} max={s['tau_rebound_max']}")
        print(f"        mech P1/P_multi/P2/no = {s['frac_P1']}/{s['frac_P_multi']}/"
              f"{s['frac_P2_new']}/{s['frac_no_switch']}  "
              f"gap-closing={s['gap_closing_share_contraction']} "
              f"(d_act={s['d_active_median']}, d_comp={s['d_competing_median']})")

    # ---- alignment CSV ----
    al_rows = []
    for K in K_list:
        for s, vals in sorted(res[K]["align"].items()):
            al_rows.append({
                "K": K, "s": s, "n": len(vals),
                "B_median": f"{float(np.median(vals)):+.6f}",
                "B_p25": f"{float(np.percentile(vals, 25)):+.6f}",
                "B_p75": f"{float(np.percentile(vals, 75)):+.6f}",
                "share_gate_open": f"{float(np.mean(vals > 0)):.4f}",
            })
    write_csv(out_dir / "gate_alignment.csv", list(al_rows[0]), al_rows)
    print("\n=== B aligned at rebound (K=100) ===")
    for r in al_rows:
        if r["K"] == 100:
            print(f"  t_r{int(r['s']):+d}: B_med={r['B_median']}  "
                  f"gate_open={r['share_gate_open']}")

    # ---- terminal vs ordinary ----
    term_rows = []
    for K in K_list:
        eps = res[K]["episodes"]
        term = next(e for e in eps if e["terminal"])
        ord_len = [e["tau_rebound"] for e in eps if not e["terminal"]]
        ord_ret = [e["tau_return"] for e in eps if not e["terminal"] and e["tau_return"] > 0]
        term_rows.append({
            "K": K, "terminal_len": term["tau_rebound"],
            "terminal_tau_return": term["tau_return"],
            "terminal_B_pos_rounds":
                int(sum(1 for s in range(1, term["tau_rebound"])
                        if float(res[K]["B"][term["t_r"] + s]) > 0)),
            "ordinary_len_median": f"{np.median(ord_len):.2f}",
            "ordinary_tau_return_median": f"{np.median(ord_ret):.2f}",
            "note": "terminal episode: last rebound -> exact stop (no next rebound)",
        })
    write_csv(out_dir / "terminal_comparison.csv", list(term_rows[0]), term_rows)
    print("\n=== terminal vs ordinary episodes ===")
    for r in term_rows:
        print(f"  K={r['K']:>3}: terminal len={r['terminal_len']} "
              f"tau_return={r['terminal_tau_return']} "
              f"B>0 rounds in episode={r['terminal_B_pos_rounds']} | "
              f"ordinary len med={r['ordinary_len_median']} "
              f"tau_return med={r['ordinary_tau_return_median']}")

    # ---- figures ----
    plot_gate_episode_alignment(res, K_list, path=out_dir / "fig_gate_episode_alignment.png")
    plot_boundary_map(res, K_primary=100, path=out_dir / "fig_boundary_map_K100.png")
    plot_peak_gap_alignment(res, K_primary=100, path=out_dir / "fig_peak_gap_alignment_K100.png")
    plot_terminal_vs_ordinary(res, K_list, path=out_dir / "fig_terminal_vs_ordinary.png")

    # ---- metadata ----
    output_files = sorted(p.name for p in out_dir.iterdir()
                          if p.is_file() and p.name != "metadata.json")
    metadata = {
        "model": cfg["model"], "model_version": cfg["model_version"],
        "experiment": cfg["experiment"],
        "run_date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": cfg,
        "headline": {
            "worst_pair_update_identity_residual":
                max(res[K]["worst_pair_resid"] for K in K_list),
            "gate_never_returns_share": {str(K): s["return_never"]
                                         for K, s in zip(K_list, sum_rows)},
            "gap_closing_share_contraction": {str(K): s["gap_closing_share_contraction"]
                                              for K, s in zip(K_list, sum_rows)},
        },
        "environment": {"python": sys.version, "platform": platform.platform(),
                        "numpy": np.__version__},
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m1c6_boundary_return/run_analysis.py",
        "outputs": {name: sha256_of(out_dir / name) for name in output_files},
    }
    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\noutputs written to: {out_dir} ({len(output_files) + 1} files)")
    print(f"runtime: {metadata['runtime_seconds']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
