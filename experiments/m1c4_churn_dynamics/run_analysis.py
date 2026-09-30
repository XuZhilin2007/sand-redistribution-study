"""Runner for M1C.4: churn-regime scaling & switching audit.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m1c4_churn_dynamics/run_analysis.py

AUDIT of the existing canonical Variant A deterministic trajectories
(K in {50,100,200,400}, q_near, alpha=0.25, exact stopping). No model
change, no new mechanism experiment.

Per active round the runner computes the M1C.3 rebound-pressure ratio

    R_t = max_{j: gamma_t(j) > 0} M_t max(Delta(j), 0) / gamma_t(j),

which satisfies R_t > 1 iff the round's D_max rebounds (exact criterion,
M1C.3); equivalence is re-verified on every round. Analyses:

  * churn window (M1B.2 operational divider: first D_max < 0.05 through
    t_stop - 1) with cross-K quantiles of R, a_t, D_max, moved mass,
    |Delta a|; rebound fractions; run lengths and empirical transition
    frequencies of the rebound/contraction symbol sequence;
  * symbol lag agreement + autocorrelation (R_t, D_max) up to lag 100;
  * low-dimensional near-return audit on standardized summaries
    (D_max, a_t, M_t, R_t), stride-sampled;
  * boundary-switch audit (|Delta a| > 0.2, M1A.1 convention) vs rebounds;
  * terminal audit (last 5/10/20 rounds vs churn bulk) and cross-K
    terminal alignment over tau = t_stop - t.
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
from sand_m0.kernel_theory import target_matched_excess  # noqa: E402
from sand_m0.model import DISTRIBUTIONS  # noqa: E402
from sand_m0.plotting import (  # noqa: E402
    plot_R_threshold,
    plot_R_crossK,
    plot_phase_a_vs_R,
    plot_terminal_alignment,
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


def fmt(v, spec: str = ".10e") -> str:
    return format(v, spec) if isinstance(v, float) else v


def analyze_run(K: int, run: dict, alpha: float, g_bins: np.ndarray,
                cfg: dict) -> dict:
    """Per-round churn diagnostics + aggregate summaries for one K."""
    stop = run["stopped_at"]
    rows = run["rows"]
    states = np.asarray(run["states"], dtype=float)
    n = stop
    T = np.arange(1, K + 1) / K
    delta = np.cumsum(g_bins) - T          # Delta(j) = G(j) - T(j), uniform target

    per = []
    mismatch_R = 0
    raw = {"R": [], "a": [], "M": [], "D": [], "da": [], "sw": []}
    for t in range(n):
        mass = states[t]
        D = cumulative_excess(mass, K)
        d_max = float(D.max())
        j_star = int(np.argmax(D)) + 1
        moved = alpha * float(mass[:j_star].sum())
        gamma = d_max - target_matched_excess(mass, j_star, alpha)
        pos_reinj = moved * np.maximum(delta, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(gamma > 0.0, pos_reinj / np.where(gamma > 0.0, gamma, 1.0), 0.0)
        R = float(ratio.max())
        j_drive = int(np.argmax(ratio)) + 1
        rebound_excess = bool((pos_reinj - gamma).max() > 0.0)   # M1C.3 criterion
        rebound_R = R > 1.0
        if rebound_excess != rebound_R:
            mismatch_R += 1
        a_t = j_star / K
        if t > 0:
            j_prev = rows[t - 1]["j_t"]
            dj = j_star - j_prev
            da = dj / K
        else:
            dj = da = None
        raw["R"].append(R)
        raw["a"].append(a_t)
        raw["M"].append(moved)
        raw["D"].append(d_max)
        raw["da"].append(abs(da) if da is not None else 0.0)
        raw["sw"].append(bool(abs(da) > cfg["switch_convention"]["threshold"])
                         if da is not None else False)
        per.append({
            "K": K, "t": t,
            "t_over_K": fmt(t / K, ".5f"),
            "t_over_tstop": fmt(t / stop, ".6f"),
            "D_max": fmt(d_max),
            "U_density": fmt(K * rows[t]["U_L2"]),
            "j_star": j_star, "a_t": fmt(a_t, ".5f"),
            "F_a": fmt(float(mass[:j_star].sum()), ".8f"),
            "moved_mass": fmt(moved, ".8f"),
            "j_drive": j_drive, "a_drive": fmt(j_drive / K, ".5f"),
            "delta_drive": fmt(float(delta[j_drive - 1]), ".8f"),
            "gamma_drive": fmt(float(gamma[j_drive - 1]), ".8f"),
            "gamma_min": fmt(float(gamma.min()), ".8f"),
            "R_t": fmt(R),
            "rebound": rebound_excess,
            "next_rebound": "",
            "j_next": rows[t + 1]["j_t"] if t + 1 < len(rows) else "",
            "jump_dj": dj, "jump_da": fmt(da, ".5f") if da is not None else "",
            "large_switch": (abs(da) > cfg["switch_convention"]["threshold"]) if da is not None else "",
            "tau_to_stop": stop - t,
        })
    for t in range(n - 1):
        per[t]["next_rebound"] = per[t + 1]["rebound"]

    # churn window: first D_max < onset threshold .. last active round
    Dm = np.array(raw["D"])
    onset = int(np.nonzero(Dm < cfg["churn_window_definition"]["onset_below"])[0][0])
    win = per[onset:n]
    R_w = np.array(raw["R"][onset:])
    a_w = np.array(raw["a"][onset:])
    M_w = np.array(raw["M"][onset:])
    da_w = np.array(raw["da"][onset:])
    reb_w = np.array([r["rebound"] for r in win], dtype=bool)

    def q(vals: np.ndarray) -> dict:
        return {f"p{p}": float(np.percentile(vals, p)) for p in [5, 25, 50, 75, 95]}

    qR, qA, qM, qDA = q(R_w), q(a_w), q(M_w), q(da_w)
    # run lengths
    runs_reb, runs_con, cur = [], [], 1
    for i in range(1, len(reb_w)):
        if reb_w[i] == reb_w[i - 1]:
            cur += 1
        else:
            (runs_reb if reb_w[i - 1] else runs_con).append(cur)
            cur = 1
    (runs_reb if reb_w[-1] else runs_con).append(cur)
    # empirical transition frequencies (NOT Markov probabilities)
    trans = {k: 0 for k in ["R->R", "R->C", "C->R", "C->C"]}
    for i in range(1, len(reb_w)):
        if reb_w[i] and reb_w[i - 1]:
            trans["R->R"] += 1
        elif not reb_w[i] and not reb_w[i - 1]:
            trans["C->C"] += 1
        elif reb_w[i - 1]:
            trans["R->C"] += 1
        else:
            trans["C->R"] += 1

    # symbol lag agreement + autocorrelations
    lag_rows = []
    s = reb_w.astype(float)
    R_c = R_w - R_w.mean()
    D_w = Dm[onset:]
    D_c = D_w - D_w.mean()
    denom_R = float(np.dot(R_c, R_c)) or 1.0
    denom_D = float(np.dot(D_c, D_c)) or 1.0
    for lag in range(1, cfg["symbol_lag_analysis"]["max_lag"] + 1):
        agree = float((s[lag:] == s[:-lag]).mean())
        acR = float(np.dot(R_c[lag:], R_c[:-lag]) / denom_R)
        acD = float(np.dot(D_c[lag:], D_c[:-lag]) / denom_D)
        lag_rows.append({"K": K, "lag": lag, "symbol_agreement": fmt(agree, ".6f"),
                         "R_autocorr": fmt(acR, ".6f"), "D_max_autocorr": fmt(acD, ".6f")})

    # low-dim near-return audit (stride-sampled summaries)
    S = np.column_stack([Dm[onset:], a_w, M_w, R_w])
    sd = S.std(axis=0)
    sd[sd == 0] = 1.0
    S = (S - S.mean(axis=0)) / sd
    stride = max(1, S.shape[0] // cfg["recurrence_audit"]["stride_cap"])
    S_s = S[::stride]
    excl = cfg["recurrence_audit"]["sampled_lag_exclusion"]
    best = (np.inf, -1, -1)
    dists = []
    for i in range(excl, len(S_s)):
        d = np.sqrt(((S_s[:i - excl] - S_s[i]) ** 2).sum(axis=1))
        dmin = float(d.min()) if len(d) else np.inf
        dists.append(dmin)
        if dmin < best[0]:
            best = (dmin, i * stride + onset, int(np.argmin(d)) * stride + onset)
    dists = np.array(dists) if dists else np.array([np.inf])
    recurrence = {
        "K": K, "window_rounds": len(win), "stride": stride,
        "sampled_points": len(S_s),
        "nearest_return_min": fmt(float(dists.min()), ".6f"),
        "nearest_return_median": fmt(float(np.median(dists)), ".6f"),
        "example_t": best[1], "example_t_prev": best[2],
        "note": "standardized 4-dim summary; distances are upper bounds on "
                "full-density near-returns (stride sampling)",
    }

    # boundary switching audit
    sw = np.array(raw["sw"][onset:])
    reb_d = da_w[reb_w]
    con_d = da_w[~reb_w]
    p_switch_after_reb = float(sw[1:][reb_w[:-1]].mean()) if reb_w[:-1].any() else float("nan")
    p_switch_base = float(sw.mean())

    summary = {
        "K": K, "t_stop": stop, "t_onset_0.05": onset,
        "onset_over_K": fmt(onset / K, ".4f"),
        "window_rounds": len(win),
        "rebound_count": int(reb_w.sum()),
        "rebound_fraction": fmt(float(reb_w.mean()), ".5f"),
        "R_p5": fmt(qR["p5"], ".5f"), "R_p25": fmt(qR["p25"], ".5f"),
        "R_p50": fmt(qR["p50"], ".5f"), "R_p75": fmt(qR["p75"], ".5f"),
        "R_p95": fmt(qR["p95"], ".5f"),
        "a_p5": fmt(qA["p5"], ".5f"), "a_p50": fmt(qA["p50"], ".5f"),
        "a_p95": fmt(qA["p95"], ".5f"),
        "D_max_p5": fmt(float(np.percentile(D_w, 5)), ".6f"),
        "D_max_p50": fmt(float(np.percentile(D_w, 50)), ".6f"),
        "D_max_p95": fmt(float(np.percentile(D_w, 95)), ".6f"),
        "M_p5": fmt(qM["p5"], ".5f"), "M_p50": fmt(qM["p50"], ".5f"),
        "M_p95": fmt(qM["p95"], ".5f"),
        "absda_p50": fmt(qDA["p50"], ".5f"), "absda_p95": fmt(qDA["p95"], ".5f"),
        "switch_freq": fmt(p_switch_base, ".5f"),
        "switch_after_rebound_freq": fmt(p_switch_after_reb, ".5f"),
        "run_reb_median": fmt(float(np.median(runs_reb)), ".3f"),
        "run_reb_max": max(runs_reb),
        "run_con_median": fmt(float(np.median(runs_con)), ".3f"),
        "run_con_max": max(runs_con),
        "RR": trans["R->R"], "RC": trans["R->C"],
        "CR": trans["C->R"], "CC": trans["C->C"],
        "criterion_mismatches": mismatch_R,
    }
    reb_da_detail = {
        "K": K,
        "absda_reb_p50": fmt(float(np.median(reb_d)), ".5f"),
        "absda_con_p50": fmt(float(np.median(con_d)), ".5f"),
        "absda_reb_p95": fmt(float(np.percentile(reb_d, 95)), ".5f"),
        "absda_con_p95": fmt(float(np.percentile(con_d, 95)), ".5f"),
    }
    return {"per": per, "summary": summary, "lag_rows": lag_rows,
            "recurrence": recurrence, "switch_detail": reb_da_detail,
            "onset": onset, "stop": stop, "R_w": R_w, "a_w": a_w,
            "reb_w": reb_w, "Dm_w": D_w, "M_w": M_w}


def terminal_audit(res: dict, cfg: dict) -> tuple[list[dict], list[dict]]:
    """Last-5/10/20 rounds vs churn bulk + cross-K tau alignment."""
    term_rows, align_rows = [], []
    for K in [50, 100, 200, 400]:
        per = [r for r in res[(K, "per")]]
        bulk_R = res[(K, "R_w")]
        bulk_reb = res[(K, "reb_w")]
        for w in cfg["terminal_windows_rounds"]:
            seg = per[-w:]
            term_rows.append({
                "K": K, "window": f"last{w}",
                "R_median": fmt(float(np.median([float(r["R_t"]) for r in seg])), ".5f"),
                "rebound_frac": fmt(float(np.mean([r["rebound"] for r in seg])), ".5f"),
                "M_median": fmt(float(np.median([float(r["moved_mass"]) for r in seg])), ".5f"),
                "D_max_median": fmt(float(np.median([float(r["D_max"]) for r in seg])), ".6f"),
                "a_median": fmt(float(np.median([float(r["a_t"]) for r in seg])), ".4f"),
                "a_min": fmt(float(np.min([float(r["a_t"]) for r in seg])), ".4f"),
                "a_max": fmt(float(np.max([float(r["a_t"]) for r in seg])), ".4f"),
            })
        term_rows.append({
            "K": K, "window": "churn_bulk",
            "R_median": fmt(float(np.median(bulk_R)), ".5f"),
            "rebound_frac": fmt(float(bulk_reb.mean()), ".5f"),
            "M_median": fmt(float(np.median(res[(K, "M_w")])), ".5f"),
            "D_max_median": fmt(float(np.median(res[(K, "Dm_w")])), ".6f"),
            "a_median": fmt(float(np.median(res[(K, "a_w")])), ".4f"),
            "a_min": fmt(float(np.min(res[(K, "a_w")])), ".4f"),
            "a_max": fmt(float(np.max(res[(K, "a_w")])), ".4f"),
        })
        for tau in range(1, cfg["terminal_alignment_tau"] + 1):
            t = res[(K, "stop")] - tau
            if t < 0:
                continue
            r = per[t]
            align_rows.append({
                "K": K, "tau": tau, "t": t,
                "D_max": r["D_max"], "R_t": r["R_t"], "a_t": r["a_t"],
                "moved_mass": r["moved_mass"], "rebound": r["rebound"],
            })
    return term_rows, align_rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=Path(__file__).resolve().parent / "config.json")
    args = parser.parse_args(argv)
    config_path = args.config.resolve()
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    out_dir = config_path.parent / cfg["output_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    K_list = [50, 100, 200, 400]
    alpha, tol = 0.25, 1e-12
    q_near = DISTRIBUTIONS["near"]
    started = time.time()

    print(f"M1C.4 churn-regime audit — config: {config_path}")
    expected_stops = {50: 324, 100: 699, 200: 1445, 400: 2876}

    res: dict = {}
    for K in K_list:
        run = run_m1a_deterministic(q_near, K=K, T_max=5000, alpha=alpha,
                                    stop_tolerance=tol)
        assert run["stopped_at"] == expected_stops[K], "baseline drift"
        one = analyze_run(K, run, alpha, q_near.bin_probs(K), cfg)
        for key in ["per", "summary", "lag_rows", "recurrence",
                    "switch_detail", "onset", "stop", "R_w", "a_w",
                    "reb_w", "Dm_w", "M_w"]:
            res[(K, key)] = one[key]
        s = res[(K, "summary")]
        print(f"  [K={K}] stop={s['t_stop']} onset(0.05)={s['t_onset_0.05']} "
              f"window={s['window_rounds']} rebounds={s['rebound_count']} "
              f"({float(s['rebound_fraction']):.4f}) mismatches={s['criterion_mismatches']}")

    # ---- per-round CSV ----
    per_rows = [r for K in K_list for r in res[(K, "per")]]
    write_csv(out_dir / "per_round_churn_diagnostics.csv",
              list(per_rows[0]), per_rows)

    # ---- churn-window summary (cross-K) ----
    sum_rows = [res[(K, "summary")] for K in K_list]
    write_csv(out_dir / "churn_window_summary.csv", list(sum_rows[0]), sum_rows)
    print("\n=== churn-window summary (cross-K) ===")
    for s in sum_rows:
        print(f"  K={s['K']:>3}: window={s['window_rounds']:>4} reb={float(s['rebound_fraction']):.4f} "
              f"R p5/50/95 = {float(s['R_p5']):.3f}/{float(s['R_p50']):.3f}/{float(s['R_p95']):.3f}  "
              f"a p50 = {float(s['a_p50']):.3f}  M p50 = {float(s['M_p50']):.4f}  "
              f"|da| p50 = {float(s['absda_p50']):.3f}")
    print("\n=== switching (empirical frequencies, deterministic) ===")
    for s in sum_rows:
        print(f"  K={s['K']:>3}: runs R med/max = {s['run_reb_median']}/{s['run_reb_max']}  "
              f"C med/max = {s['run_con_median']}/{s['run_con_max']}  "
              f"RR/RC/CR/CC = {s['RR']}/{s['RC']}/{s['CR']}/{s['CC']}  "
              f"switch_freq={float(s['switch_freq']):.3f} "
              f"P(switch|rebound prev)={float(s['switch_after_rebound_freq']):.3f}")

    # ---- switching detail / lag agreement / recurrence ----
    sw_rows = [res[(K, "switch_detail")] for K in K_list]
    write_csv(out_dir / "boundary_switch_summary.csv", list(sw_rows[0]), sw_rows)
    lag_rows = [r for K in K_list for r in res[(K, "lag_rows")]]
    write_csv(out_dir / "lag_agreement.csv", list(lag_rows[0]), lag_rows)
    rec_rows = [res[(K, "recurrence")] for K in K_list]
    write_csv(out_dir / "recurrence_summary.csv", list(rec_rows[0]), rec_rows)
    print("\n=== periodicity / recurrence ===")
    for K, r in zip(K_list, rec_rows):
        lags = [x for x in res[(K, "lag_rows")]]
        best_lag = max(lags, key=lambda x: float(x["symbol_agreement"]))
        print(f"  K={K:>3}: nearest-return min/median = {r['nearest_return_min']}/"
              f"{r['nearest_return_median']} (stride {r['stride']}); "
              f"max symbol agreement {float(best_lag['symbol_agreement']):.3f} "
              f"at lag {best_lag['lag']}")

    # ---- terminal audit + alignment ----
    term_rows, align_rows = terminal_audit(res, cfg)
    write_csv(out_dir / "terminal_summary.csv", list(term_rows[0]), term_rows)
    write_csv(out_dir / "terminal_alignment.csv", list(align_rows[0]), align_rows)
    print("\n=== terminal audit (last 5/10/20 vs churn bulk) ===")
    for K in K_list:
        rs = [r for r in term_rows if r["K"] == K]
        for r in rs:
            print(f"  K={K:>3} {r['window']:>10}: R_med={r['R_median']} "
                  f"reb={r['rebound_frac']} M_med={r['M_median']} "
                  f"D_max_med={r['D_max_median']} a range [{r['a_min']},{r['a_max']}]")

    # ---- figures ----
    plot_R_threshold(res, K_primary=100, path=out_dir / "fig_R_threshold_K100.png")
    plot_R_crossK(res, K_list, path=out_dir / "fig_R_crossK_quantiles.png")
    plot_phase_a_vs_R(res, K_primary=100, path=out_dir / "fig_phase_a_vs_R_K100.png")
    plot_terminal_alignment(align_rows, path=out_dir / "fig_terminal_alignment.png")

    # ---- metadata ----
    output_files = sorted(p.name for p in out_dir.iterdir()
                          if p.is_file() and p.name != "metadata.json")
    metadata = {
        "model": cfg["model"], "model_version": cfg["model_version"],
        "experiment": cfg["experiment"],
        "run_date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": cfg,
        "headline": {
            "criterion_mismatches_total":
                sum(res[(K, "summary")]["criterion_mismatches"] for K in K_list),
            "rebound_fractions_churn_window": {
                str(K): res[(K, "summary")]["rebound_fraction"] for K in K_list},
            "baseline_stops": expected_stops,
        },
        "environment": {"python": sys.version, "platform": platform.platform(),
                        "numpy": np.__version__},
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m1c4_churn_dynamics/run_analysis.py",
        "outputs": {name: sha256_of(out_dir / name) for name in output_files},
    }
    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\noutputs written to: {out_dir} ({len(output_files) + 1} files)")
    print(f"runtime: {metadata['runtime_seconds']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
