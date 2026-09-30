"""Runner for M2 Gate 2: sign-placement generality experiment.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m2_gate2_sign_placement/run_experiment.py

Research question (frozen): keeping the absolute sign-changing mismatch
profile identical (Delta_O = -Delta_A pointwise), does placing the
positive mismatch over the historical canonical witness-support region
(XA, Witness-Aligned) versus outside it (XO, Witness-Opposed) change
rebound / recurrence / churn-like dynamics on the canonical M1 dynamics?

Frozen kernels (config.json):

  * XA -- Delta_A = x/2 | 1/3 - x/2 | 2/3 - x | x - 1  (positive on
          (0, 2/3) over the canonical witness band, negative in Far)
  * XO -- Delta_O = -Delta_A pointwise (exact sign reflection; negative
          on the witness band, positive in Far)

Both share breakpoints {1/3, 2/3, 5/6}, piecewise complexity, max
|Delta| = 1/6, the canonical initial state, the canonical controller and
exact stopping. The ONLY manipulated variable is sign placement.

Fixed resolution matrix K in {50, 100, 200, 400} (no sequential
decisions). Horizon H(K) = 10*K active rounds or exact stop; censored
observations are recorded as such and never extended.

Per-round diagnostics include the Gate-1-audit witness/risk machinery
from round one: E_t(j) = M*Delta - gamma via the repository
`rebound_residual` helper, raw and prefix-restricted argmax (smallest-
index ties), sign of Delta at the risk point, and positive/negative
support classification for the migration analysis.
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

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from sand_m0.adaptive import run_m1a_deterministic  # noqa: E402
from sand_m0.kernel_theory import (  # noqa: E402
    rebound_residual,
    target_matched_excess,
)
from sand_m0.model import (  # noqa: E402
    DISTRIBUTIONS,
    witness_aligned_distribution,
    witness_opposed_distribution,
)

KERNEL_ORDER = ["XA", "XO"]
K_LIST = [50, 100, 200, 400]
COLORS = {"XA": "#1f77b4", "XO": "#d62728", "KC": "#777777"}


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_csv_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


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


def region_of(x: float) -> str:
    if x < 1.0 / 3.0:
        return "Near"
    if x <= 2.0 / 3.0:
        return "Middle"
    return "Far"


def witness_argmax(E: np.ndarray) -> tuple[int, float]:
    """Raw witness/risk argmax: smallest index on ties (np.argmax)."""
    j = int(np.argmax(E))
    return j + 1, float(E[j])


def witness_argmax_prefix(E: np.ndarray, j_star: int) -> tuple[int, float]:
    """Prefix-restricted argmax (j <= j*): criterion-relevant supplement
    (raw argmax can land on structurally-impossible suffix positions)."""
    j = int(np.argmax(E[:j_star]))
    return j + 1, float(E[j])


def sign_class(value: float, tau: float) -> str:
    if value > tau:
        return "positive"
    if value < -tau:
        return "negative"
    return "zero"


BRANCHES = {
    "XA": [lambda v: 1.5 * v,
           lambda v: 0.5 * v + 1.0 / 3.0,
           lambda v: 2.0 / 3.0 + 0.0 * v,
           lambda v: 2.0 * v - 1.0],
    "XO": [lambda v: 0.5 * v,
           lambda v: 1.5 * v - 1.0 / 3.0,
           lambda v: 2.0 * v - 2.0 / 3.0,
           lambda v: 1.0 + 0.0 * v],
}
BREAKPOINTS = [(1.0 / 3.0, 0, 1, "1/3"), (2.0 / 3.0, 1, 2, "2/3"),
               (5.0 / 6.0, 2, 3, "5/6")]


def static_kernel_checks(name: str, respray, other, K: int, tau: float
                         ) -> list[dict]:
    """Frozen static checks incl. the critical sign-reflection property."""
    edges = np.linspace(0.0, 1.0, K + 1)
    G_edges = respray.cdf(edges)
    raw_q = np.diff(G_edges)
    x = np.arange(1, K + 1) / K
    delta = respray.cdf(x) - x
    grid = np.linspace(0.0, 1.0, 12001)
    G_grid = respray.cdf(grid)

    def row(check: str, ok: bool, detail: str) -> dict:
        return {"kernel": name, "K": K, "check": check, "detail": detail,
                "pass": bool(ok)}

    rows = [
        row("G(0)=0", abs(float(respray.cdf(np.array([0.0]))[0])) <= tau,
            f"G(0)={float(respray.cdf(np.array([0.0]))[0]):.3e}"),
        row("G(1)=1", abs(float(respray.cdf(np.array([1.0]))[0]) - 1.0) <= tau,
            f"G(1)={float(respray.cdf(np.array([1.0]))[0]):.16f}"),
        row("G nondecreasing (bin edges)", bool(np.all(np.diff(G_edges) >= -tau)),
            f"min diff={float(np.diff(G_edges).min()):.3e}"),
        row("G nondecreasing (dense grid)", bool(np.all(np.diff(G_grid) >= -tau)),
            f"min diff={float(np.diff(G_grid).min()):.3e}"),
        row("q_j >= -tau_num", float(raw_q.min()) >= -tau,
            f"min q={float(raw_q.min()):.3e}"),
        row("|sum(q)-1| <= tau_num", abs(float(raw_q.sum()) - 1.0) <= tau,
            f"sum={float(raw_q.sum()):.16e}"),
        row("density pdf >= 0 on dense grid",
            bool(np.all(respray.pdf(grid) >= -tau)),
            f"min pdf={float(respray.pdf(grid).min()):.3e}"),
        row("sign reflection |Delta_O + Delta_A| <= tau on grid",
            float(np.abs(delta + (other.cdf(x) - x)).max()) <= tau,
            f"max={float(np.abs(delta + (other.cdf(x) - x)).max()):.3e}"),
        row("sign reflection G_O + G_A == 2x on grid",
            float(np.abs(respray.cdf(x) + other.cdf(x) - 2.0 * x).max()) <= tau,
            f"max={float(np.abs(respray.cdf(x) + other.cdf(x) - 2.0 * x).max()):.3e}"),
        row("extrema: |Delta(1/3)| = |Delta(5/6)| = 1/6 with frozen signs "
            "(XA: +1/6 @1/3, -1/6 @5/6; XO reflected); grid bounded by them",
            abs(abs(float(respray.cdf(np.array([1.0 / 3.0]))[0]) - 1.0 / 3.0)
                - 1.0 / 6.0) <= tau
            and abs(abs(float(respray.cdf(np.array([5.0 / 6.0]))[0]) - 5.0 / 6.0)
                    - 1.0 / 6.0) <= tau
            and sign_class(float(respray.cdf(np.array([1.0 / 3.0]))[0]) - 1.0 / 3.0,
                           tau) == ("positive" if name == "XA" else "negative")
            and sign_class(float(respray.cdf(np.array([5.0 / 6.0]))[0]) - 5.0 / 6.0,
                           tau) == ("negative" if name == "XA" else "positive")
            and float(delta.max()) <= 1.0 / 6.0 + tau
            and float(delta.min()) >= -1.0 / 6.0 - tau,
            f"Delta(1/3)={float(respray.cdf(np.array([1.0 / 3.0]))[0]) - 1.0 / 3.0:.10f}, "
            f"Delta(5/6)={float(respray.cdf(np.array([5.0 / 6.0]))[0]) - 5.0 / 6.0:.10f}, "
            f"grid range=[{float(delta.min()):.6f},{float(delta.max()):.6f}]"),
        row("sign supports (Delta>0 on own positive support)",
            all(sign_class(float(respray.cdf(np.array([v]))[0]) - v, tau)
                == s for v, s in [(0.2, "positive" if name == "XA" else "negative"),
                                  (0.8, "negative" if name == "XA" else "positive")]),
            "probe x=0.2, x=0.8"),
    ]
    for bp, i_l, i_r, tag in BREAKPOINTS:
        # branch formulas evaluated AT the breakpoint (a bp±h probe would
        # measure O(h) slope difference, not continuity)
        g_fun = float(respray.cdf(np.array([bp]))[0])
        g_l = float(BRANCHES[name][i_l](bp))
        g_r = float(BRANCHES[name][i_r](bp))
        gap = max(abs(g_fun - g_l), abs(g_fun - g_r), abs(g_l - g_r))
        rows.append(row(f"G continuous at {tag}", gap <= tau,
                        f"branch gap={gap:.3e}"))
    return rows


def per_round_diagnostics(kernel: str, K: int, run: dict, g_bins: np.ndarray,
                          respray, alpha: float, tau: float,
                          switch_threshold: float) -> tuple[list[dict], dict]:
    """Per-round diagnostics with witness/risk + sign + support fields
    (E via the repository rebound_residual helper)."""
    rows = run["rows"]
    stop = run["stopped_at"]
    n = stop if stop is not None else len(rows) - 1
    states = np.asarray(run["states"], dtype=float)
    T = np.arange(1, K + 1) / K
    delta = np.cumsum(g_bins) - T

    out: list[dict] = []
    tp = tn = fp = fn = n_boundary = fp_nb = fn_nb = 0
    max_resid = max_mass_err = 0.0
    min_bin = np.inf
    for t in range(n):
        mass = states[t]
        row = rows[t]
        j_star = int(row["j_t"])
        d0 = float(row["D_max"])
        moved = float(row["removed_mass"])
        D_tm = target_matched_excess(mass, j_star, alpha)
        gamma = d0 - D_tm
        E = rebound_residual(mass, j_star, alpha, g_bins)
        j_risk, worst = witness_argmax(E)
        j_pref, worst_pref = witness_argmax_prefix(E, j_star)
        x_risk, x_pref = j_risk / K, j_pref / K
        pos_reinj = moved * np.maximum(delta, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(gamma > 0.0,
                             pos_reinj / np.where(gamma > 0.0, gamma, 1.0), 0.0)
        R = float(ratio.max())
        criterion_rebound = worst > 0.0
        pred_profile = D_tm + moved * delta
        actual_profile = np.cumsum(states[t + 1]) - T
        resid = float(np.abs(pred_profile - actual_profile).max())
        max_resid = max(max_resid, resid)
        d_after = float(rows[t + 1]["D_max"])
        actual_rebound = d_after > d0
        boundary_case = abs(worst) <= tau
        n_boundary += int(boundary_case)
        if actual_rebound and criterion_rebound:
            tp += 1
        elif not actual_rebound and not criterion_rebound:
            tn += 1
        elif criterion_rebound:
            fp += 1
            fp_nb += int(not boundary_case)
        else:
            fn += 1
            fn_nb += int(not boundary_case)
        max_mass_err = max(max_mass_err, abs(float(states[t + 1].sum()) - 1.0))
        min_bin = min(min_bin, float(states[t + 1].min()))
        da = (j_star - rows[t - 1]["j_t"]) / K if t > 0 else None
        delta_at_risk = float(respray.cdf(np.array([x_risk]))[0]) - x_risk
        delta_at_pref = float(respray.cdf(np.array([x_pref]))[0]) - x_pref
        out.append({
            "kernel": kernel, "K": K, "t": t,
            "D_max_before": f"{d0:.10e}",
            "j_t": j_star, "a_t": f"{j_star / K:.5f}",
            "moved_mass": f"{moved:.10e}",
            "D_max_after": f"{d_after:.10e}",
            "d_D_max": f"{d_after - d0:.10e}",
            "rebound": bool(actual_rebound),
            "delta_a": f"{da:.5f}" if da is not None else "",
            "large_switch": (abs(da) > switch_threshold) if da is not None else "",
            "decomposition_residual": f"{resid:.6e}",
            "criterion_rebound": bool(criterion_rebound),
            "worst_excess": f"{worst:.8e}",
            "numerical_boundary": bool(boundary_case),
            "R_t": f"{R:.10e}",
            "j_risk": j_risk,
            "x_risk": f"{x_risk:.5f}",
            "argmax_on_suffix": bool(j_risk > j_star),
            "j_risk_prefix": j_pref,
            "x_risk_prefix": f"{x_pref:.5f}",
            "sign_at_risk": sign_class(delta_at_risk, tau),
            "sign_at_risk_prefix": sign_class(delta_at_pref, tau),
            "region_risk": region_of(x_risk),
            "region_risk_prefix": region_of(x_pref),
            "stopped": stop is not None,
            "t_stop": stop if stop is not None else "",
        })
    return out, {
        "active_rounds": n, "TP": tp, "TN": tn, "FP": fp, "FN": fn,
        "numerical_boundary": n_boundary, "FP_nonboundary": fp_nb,
        "FN_nonboundary": fn_nb, "max_decomposition_residual": max_resid,
        "max_mass_error": max_mass_err, "min_bin_mass": min_bin,
    }


def classify(N_R: int, S_R_over_K: float | None) -> str:
    if N_R == 0:
        return "No rebound"
    if N_R == 1:
        return "Isolated rebound"
    if S_R_over_K is not None and S_R_over_K >= 1.0:
        return "Long-span recurrent / churn-like transient"
    return "Short recurrent transient"


def trajectory_summary(kernel: str, K: int, run: dict, per_rows: list[dict],
                       theory: dict) -> dict:
    n = theory["active_rounds"]
    stop = run["stopped_at"]
    reb_t = [int(r["t"]) for r in per_rows if r["rebound"]]
    n_reb = len(reb_t)
    f_R = n_reb / n if n else 0.0
    if n_reb >= 1:
        t_first, t_last = reb_t[0], reb_t[-1]
        S_R = t_last - t_first
        S_R_over_K = S_R / K
        inter = np.diff(np.array(reb_t)) if n_reb >= 2 else np.array([])
    else:
        t_first = t_last = S_R = ""
        S_R_over_K = None
        inter = np.array([])
    da = np.array([abs(float(r["delta_a"])) for r in per_rows[1:]]
                  if len(per_rows) > 1 else [])
    switches = [r["large_switch"] for r in per_rows[1:]]
    moved = np.array([float(r["moved_mass"]) for r in per_rows])
    d0 = np.array([float(r["D_max_before"]) for r in per_rows])
    R_all = np.array([float(r["R_t"]) for r in per_rows])
    longest_run = cur = 0
    for r in per_rows:
        cur = 0 if r["rebound"] else cur + 1
        longest_run = max(longest_run, cur)
    return {
        "kernel": kernel, "K": K,
        "active_rounds": n,
        "stop_status": ("exact_stop" if stop is not None
                        else "not stopped within pre-registered Gate-2 horizon"),
        "t_stop": stop if stop is not None else "",
        "N_R": n_reb,
        "f_R": f"{f_R:.5f}",
        "S_R": S_R,
        "S_R_over_K": f"{S_R_over_K:.5f}" if S_R_over_K is not None else "",
        "t_first_rebound": t_first,
        "t_last_rebound": t_last,
        "rebound_interarrival_median":
            f"{float(np.median(inter)):.3f}" if len(inter) else "",
        "longest_contraction_run": longest_run,
        "boundary_switch_count": sum(1 for s in switches if s),
        "boundary_switch_freq":
            f"{sum(1 for s in switches if s) / max(len(switches), 1):.5f}",
        "abs_delta_a_median": f"{float(np.median(da)):.5f}" if len(da) else "",
        "moved_mass_median": f"{float(np.median(moved)):.6f}",
        "moved_mass_total": f"{float(moved.sum()):.6f}",
        "D_max_initial": f"{d0[0]:.6f}" if len(d0) else "",
        "D_max_max": f"{d0.max():.6f}" if len(d0) else "",
        "D_max_final": f"{run['rows'][n]['D_max']:.6e}",
        "max_R_t": f"{float(R_all.max()):.5f}" if len(R_all) else "",
        "max_decomposition_residual": f"{theory['max_decomposition_residual']:.3e}",
        "criterion_TP": theory["TP"], "criterion_TN": theory["TN"],
        "criterion_FP": theory["FP"], "criterion_FN": theory["FN"],
        "numerical_boundary_count": theory["numerical_boundary"],
        "FP_nonboundary": theory["FP_nonboundary"],
        "FN_nonboundary": theory["FN_nonboundary"],
        "max_mass_error": f"{theory['max_mass_error']:.3e}",
        "min_bin_mass": f"{theory['min_bin_mass']:.3e}",
        "qualitative_category": classify(n_reb, S_R_over_K),
    }


def support_migration_rows(per_rows: list[dict]) -> list[dict]:
    """Table 3: risk/witness location vs the kernel's own positive support,
    for all rounds / top-25%-R_t rounds / actual rebound rounds, in raw and
    prefix-restricted views."""
    out: list[dict] = []
    if not per_rows:
        return out
    kernel, K = per_rows[0]["kernel"], per_rows[0]["K"]
    R = np.array([float(r["R_t"]) for r in per_rows])
    top25_cut = float(np.percentile(R, 75))
    for view, x_field, s_field in [("raw", "x_risk", "sign_at_risk"),
                                   ("prefix", "x_risk_prefix", "sign_at_risk_prefix")]:
        subsets = {
            "all active rounds": per_rows,
            "top 25% R_t rounds": [r for r in per_rows if float(r["R_t"]) >= top25_cut],
            "actual rebound rounds": [r for r in per_rows if r["rebound"]],
        }
        for tag, subset in subsets.items():
            if not subset:
                continue
            x = np.array([float(r[x_field]) for r in subset])
            signs = [r[s_field] for r in subset]
            n = len(subset)
            out.append({
                "kernel": kernel, "K": K, "view": view, "sample": tag, "n": n,
                "share_positive": f"{signs.count('positive') / n:.5f}",
                "share_negative": f"{signs.count('negative') / n:.5f}",
                "share_zero": f"{signs.count('zero') / n:.5f}",
                "x_median": f"{float(np.median(x)):.5f}",
                "x_p05": f"{float(np.percentile(x, 5)):.5f}",
                "x_p95": f"{float(np.percentile(x, 95)):.5f}",
            })
    return out


def _k_panel_axes():
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2))
    return fig, np.asarray(axes).reshape(-1)


def plot_fig(ax_data: dict, ylabel: str, path: Path, title: str,
             logy: bool = False, mark_two_thirds: bool = False) -> None:
    fig, axes = _k_panel_axes()
    for ax, K in zip(axes, K_LIST):
        for kernel in KERNEL_ORDER:
            arr = ax_data[(kernel, K)]
            t = np.arange(len(arr["y"]))
            if logy:
                ax.semilogy(t, np.maximum(arr["y"], 1e-16), color=COLORS[kernel],
                            lw=0.9, label=kernel)
            else:
                ax.plot(t, arr["y"], color=COLORS[kernel], lw=0.8, label=kernel)
            if arr.get("reb") is not None and arr["reb"].any():
                ax.scatter(t[arr["reb"]], np.maximum(arr["y"][arr["reb"]], 1e-16)
                           if logy else arr["y"][arr["reb"]],
                           s=10, color="black", zorder=3)
        if mark_two_thirds:
            ax.axhline(2.0 / 3.0, color="gray", ls="--", lw=0.8)
            ax.text(0.01, 2.0 / 3.0, " x=2/3 sign boundary", fontsize=7,
                    va="bottom", color="gray")
        ax.set_title(f"K={K}", fontsize=10)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("round t")
        ax.legend(fontsize=7, loc="best")
    fig.suptitle(title)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=150)
    plt.close(fig)


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
    tau = float(cfg["tau_num"])
    alpha = float(cfg["alpha"])
    tol = float(cfg["stop_tolerance_exact"])
    switch_threshold = float(cfg["switch_convention"]["threshold"])
    q_near = DISTRIBUTIONS["near"]
    dists = {"XA": witness_aligned_distribution(),
             "XO": witness_opposed_distribution()}
    started = time.time()

    print(f"M2 Gate 2 sign-placement experiment — config: {config_path}")
    print(f"  alpha={alpha}, exact tol={tol}, tau_num={tau}, horizon H(K)=10K")

    # ---------------- A4: static validation (incl. sign reflection) ----------------
    print("\n=== A4: static validation (tau_num = %.0e) ===" % tau)
    static_rows: list[dict] = []
    for name in KERNEL_ORDER:
        for K in K_LIST:
            static_rows.extend(
                static_kernel_checks(name, dists[name],
                                     dists["XO" if name == "XA" else "XA"],
                                     K, tau))
    write_csv(out_dir / "static_kernel_validation.csv",
              ["kernel", "K", "check", "detail", "pass"], static_rows)
    fails = [r for r in static_rows if not r["pass"]]
    print(f"  static checks: {len(static_rows) - len(fails)}/{len(static_rows)} passed")
    if fails:
        for r in fails:
            print(f"  [FAIL] {r['kernel']} K={r['K']}: {r['check']} ({r['detail']})")
        print("STATIC VALIDATION FAILED — no dynamics will be run.")
        return 1

    # ---------------- B: fixed resolution matrix ----------------
    print("\n=== B: XA / XO @ K = 50, 100, 200, 400 (fixed matrix) ===")
    per_round_all: list[dict] = []
    summary_rows: list[dict] = []
    migration_rows: list[dict] = []
    per_by: dict[tuple[str, int], dict] = {}
    total_tp = total_tn = total_fp = total_fn = total_nb = 0
    total_fp_nb = total_fn_nb = 0
    max_resid_all = 0.0
    for kernel in KERNEL_ORDER:
        for K in K_LIST:
            g = dists[kernel].bin_probs(K)
            run = run_m1a_deterministic(
                q_near, K=K, T_max=10 * K, alpha=alpha, stop_tolerance=tol,
                redistribution="throw", redistribution_dist=dists[kernel])
            states = np.asarray(run["states"], dtype=float)
            if not np.array_equal(states[0], q_near.bin_probs(K)):
                print(f"[FAIL] shared initial state violated: {kernel} K={K}")
                return 1
            per, theory = per_round_diagnostics(kernel, K, run, g,
                                                dists[kernel], alpha, tau,
                                                switch_threshold)
            # witness sign gate: an actual rebound witness must sit on
            # positive mismatch (prefix-restricted view); a NEGATIVE
            # classification is a structural contradiction (indexing/logic
            # gate, config §29); ZERO is only tolerable on numerical-boundary
            # rounds (|worst excess| <= tau) and is recorded as a warning
            bad = [r for r in per if r["rebound"]
                   and r["sign_at_risk_prefix"] == "negative"]
            zero_w = [r for r in per if r["rebound"]
                      and r["sign_at_risk_prefix"] == "zero"]
            if bad or [r for r in zero_w if not r["numerical_boundary"]]:
                print(f"[FAIL] rebound witness on non-positive mismatch: "
                      f"{kernel} K={K} "
                      f"neg t={[r['t'] for r in bad]}, "
                      f"non-boundary zero t={[r['t'] for r in zero_w if not r['numerical_boundary']]} "
                      "— indexing/logic gate (config §29).")
                return 1
            if zero_w:
                print(f"  [warn] {kernel} K={K}: zero-sign witnesses on "
                      f"numerical-boundary rounds t={[r['t'] for r in zero_w]}")
            summary = trajectory_summary(kernel, K, run, per, theory)
            per_round_all.extend(per)
            summary_rows.append(summary)
            migration_rows.extend(support_migration_rows(per))
            per_by[(kernel, K)] = {
                "x": np.array([float(r["a_t"]) for r in per]),
                "reb": np.array([bool(r["rebound"]) for r in per]),
                "d": np.array([float(r["D_max_before"]) for r in per]
                              + [float(run["rows"][theory["active_rounds"]]["D_max"])]),
                "R": np.array([float(r["R_t"]) for r in per]),
                "xr": np.array([float(r["x_risk_prefix"]) for r in per]),
            }
            total_tp += theory["TP"]
            total_tn += theory["TN"]
            total_fp += theory["FP"]
            total_fn += theory["FN"]
            total_nb += theory["numerical_boundary"]
            total_fp_nb += theory["FP_nonboundary"]
            total_fn_nb += theory["FN_nonboundary"]
            max_resid_all = max(max_resid_all, theory["max_decomposition_residual"])
            stop = run["stopped_at"]
            print(f"  [{kernel}] K={K}: active={theory['active_rounds']}, "
                  f"stop={'t=' + str(stop) if stop is not None else 'CENSORED at 10K'}, "
                  f"N_R={summary['N_R']}, f_R={summary['f_R']}, "
                  f"S_R/K={summary['S_R_over_K'] or '-'}, maxR={summary['max_R_t']}, "
                  f"resid={summary['max_decomposition_residual']}, "
                  f"FP={theory['FP']}, FN={theory['FN']}, "
                  f"category={summary['qualitative_category']}")

    # ---------------- C: theory totals + tables + figures ----------------
    total_active = sum(r["active_rounds"] for r in summary_rows)
    print("\n=== C: theory verification totals ===")
    print(f"  active rounds={total_active}, TP={total_tp}, TN={total_tn}, "
          f"FP={total_fp}, FN={total_fn}, numerical-boundary={total_nb}, "
          f"non-boundary FP={total_fp_nb}, non-boundary FN={total_fn_nb}, "
          f"max decomposition residual={max_resid_all:.3e}")
    if total_fp_nb != 0 or total_fn_nb != 0:
        print("UNEXPLAINED (non-boundary) CRITERION FP/FN — pausing "
              "scientific interpretation (boundary-flagged FP/FN at "
              "|worst excess| <= tau_num are recorded, semantics untouched).")
        return 1
    print("  all FP/FN are floating-boundary cases (|worst excess| <= tau_num); "
          "strict inequality semantics untouched")

    write_csv(out_dir / "per_round_diagnostics.csv",
              list(per_round_all[0]), per_round_all)
    write_csv(out_dir / "trajectory_summary.csv", list(summary_rows[0]),
              summary_rows)
    write_csv(out_dir / "support_migration.csv", list(migration_rows[0]),
              migration_rows)

    plot_fig({k: {"y": v["d"], "reb": None} for k, v in per_by.items()},
             "D_max", out_dir / "fig_gate2_D_max.png",
             "M2 Gate 2 Figure A — D_max(t): XA vs XO", logy=True)
    plot_fig({k: {"y": v["x"], "reb": None} for k, v in per_by.items()},
             "boundary a_t", out_dir / "fig_gate2_boundary.png",
             "M2 Gate 2 Figure B — sweep boundary a_t: XA vs XO")
    plot_fig({k: {"y": v["R"], "reb": v["reb"]} for k, v in per_by.items()},
             "R_t", out_dir / "fig_gate2_R_t.png",
             "M2 Gate 2 Figure C — R_t with rebound markers: XA vs XO",
             logy=True)
    plot_fig({k: {"y": v["xr"], "reb": None} for k, v in per_by.items()},
             "x_risk (prefix-restricted)", out_dir / "fig_gate2_risk_location.png",
             "M2 Gate 2 Figure D — risk/witness location x_t: XA vs XO",
             mark_two_thirds=True)

    # ---------------- contextual KC reference (committed Gate 1, no re-run) ----
    print("\n=== contextual KC reference (committed Gate 1 results, no re-run) ===")
    kc_path = (REPO_ROOT / "experiments/m2_gate1_generality/results/"
               "rebound_summary.csv")
    for r in read_csv_rows(kc_path):
        if r["kernel"] == "KC":
            print(f"  KC K={r['K']}: t_stop={r['t_stop']}, N_R={r['N_R']}, "
                  f"f_R={r['f_R']}, S_R/K={r['S_R_over_K']}, "
                  f"category={r['qualitative_category']}")

    print("\n=== Table 2 (dynamics summary) ===")
    for r in summary_rows:
        print(f"  {r['kernel']} K={r['K']:>3}: {r['stop_status'][:24]:<24} "
              f"t_stop={str(r['t_stop']):>5} N_R={r['N_R']:>4} f_R={r['f_R']} "
              f"S_R/K={str(r['S_R_over_K']):>8} maxR={r['max_R_t']} "
              f"cat={r['qualitative_category']}")
    print("\n=== Table 3 (support migration, prefix-restricted view) ===")
    for r in migration_rows:
        if r["view"] == "prefix":
            print(f"  {r['kernel']} K={r['K']:>3} [{r['sample']:<24}] "
                  f"pos={r['share_positive']} neg={r['share_negative']} "
                  f"x_med={r['x_median']} p05-p95=[{r['x_p05']},{r['x_p95']}]")

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
        "headline": {
            "theory": {"active_rounds": total_active, "TP": total_tp,
                       "TN": total_tn, "FP": total_fp, "FN": total_fn,
                       "numerical_boundary": total_nb,
                       "FP_nonboundary": total_fp_nb,
                       "FN_nonboundary": total_fn_nb,
                       "max_decomposition_residual": max_resid_all},
            "dynamics": {f"{r['kernel']}@K={r['K']}": {
                "stop_status": r["stop_status"], "t_stop": r["t_stop"],
                "N_R": r["N_R"], "f_R": r["f_R"], "S_R_over_K": r["S_R_over_K"],
                "category": r["qualitative_category"]} for r in summary_rows},
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m2_gate2_sign_placement/run_experiment.py",
        "outputs": {name: sha256_of(out_dir / name) for name in output_files},
    }
    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"\noutputs written to: {out_dir}")
    print(f"runtime: {round(time.time() - started, 2)} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
