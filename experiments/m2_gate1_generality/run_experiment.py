"""Runner for M2 Gate 1: redistribution-mismatch generality experiment.

One command regenerates every CSV, figure and the metadata file:

    python experiments/m2_gate1_generality/run_experiment.py

Research question (frozen, Stage 3 / M2 Gate 1): how do mismatch sign,
mismatch amplitude and mismatch spatial shape of a fixed redistribution
law G against the uniform target affect rebound, repeated rebound,
switching, long-span churn-like transient, and exact stopping /
finite-horizon censoring on the canonical M1 dynamics?

Frozen kernel set (config.json):

  * K0 -- G0(x) = x                 zero-mismatch target-matched reference (M1C.1 U)
  * K- -- G(x) = x^2                negative mismatch mirror (theorem-backed control)
  * KW -- G_W(x) = 5/4 x - 1/4 x^2  weak positive mismatch (lambda = 1/4, frozen)
  * KC -- G_C(x) = 2x - x^2         canonical q_near anchor (integrity only)
  * KS -- piecewise flat-top        positive spatial-shape control (peak 1/4, integral 1/6)

Every kernel shares the canonical initial state (q_near bin probabilities)
and the canonical controller (alpha = 0.25, cumulative-excess selection,
prefix removal, exact stopping at D_max <= 1e-12). Horizon: H(K) = 10*K
active rounds or exact stop, whichever first; censored observations are
recorded as such and never extended.

Execution order (frozen): A2 static kernel validation -> A3/A4 K0/KC
baseline integrity vs committed M1C.1 rows -> A5 K- acceptance gate at
K=100 -> B1/B2 KW/KS at K=100 -> C KW/KS cross-K at {50, 200, 400}.
Any failed gate stops the scientific execution before the next phase.
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
from sand_m0.kernel_theory import target_matched_excess  # noqa: E402
from sand_m0.model import (  # noqa: E402
    DISTRIBUTIONS,
    flat_top_distribution,
    weak_positive_distribution,
)

KERNEL_ORDER = ["K0", "K-", "KW", "KC", "KS"]
NEW_KERNELS = ["K-", "KW", "KS"]
CANONICAL_K_LIST = [50, 100, 200, 400]
COLORS = {"K0": "#777777", "K-": "#d62728", "KW": "#ff7f0e",
          "KC": "#1f77b4", "KS": "#2ca02c"}


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


def band(v: np.ndarray) -> dict:
    return {"p5": float(np.percentile(v, 5)), "p50": float(np.percentile(v, 50)),
            "p95": float(np.percentile(v, 95))}


def build_registry(cfg: dict) -> dict:
    """Frozen Gate 1 kernel registry: respray law + K list per kernel (A1 reuse).

    The initial state is ALWAYS the canonical q_near bin probabilities
    (passed as `dist` to the untouched canonical simulator); only the
    respray law varies. KC uses redistribution_dist=None, i.e. the exact
    canonical single-distribution code path.
    """
    q_near = DISTRIBUTIONS["near"]
    return {
        "K0": {"respray": DISTRIBUTIONS["uniform"],
               "K_list": list(cfg["kernels"]["K0"]["K_list"]),
               "role": "zero-mismatch / target-matched integrity reference (M1C.1 U reuse)"},
        "K-": {"respray": DISTRIBUTIONS["far"],
               "K_list": list(cfg["kernels"]["K-"]["K_list"]),
               "role": "negative mismatch mirror (sign control, theorem-backed audit)"},
        "KW": {"respray": weak_positive_distribution(),
               "K_list": list(cfg["kernels"]["KW"]["K_list"]),
               "role": "weak positive mismatch, amplitude ablation (lambda=1/4)"},
        "KC": {"respray": None,
               "K_list": list(cfg["kernels"]["KC"]["K_list"]),
               "role": "canonical q_near calibration anchor (integrity only)"},
        "KS": {"respray": flat_top_distribution(),
               "K_list": list(cfg["kernels"]["KS"]["K_list"]),
               "role": "positive spatial-shape control (flat-top vs parabolic hump)"},
    }


# ------------------------------- A2: static validation -------------------------------

def static_kernel_checks(name: str, respray, K: int, tau: float) -> list[dict]:
    """Frozen static checks (config §8): CDF validity + normalization +
    kernel-specific mismatch checks, all at tolerance tau_num."""
    edges = np.linspace(0.0, 1.0, K + 1)
    G_edges = respray.cdf(edges)
    raw_q = np.diff(G_edges)
    q_norm = raw_q / raw_q.sum()
    grid = np.linspace(0.0, 1.0, 10001)
    G_grid = respray.cdf(grid)
    x = np.arange(1, K + 1) / K
    delta = respray.cdf(x) - x

    def row(check: str, ok: bool, detail: str) -> dict:
        return {"kernel": name, "K": K, "check": check, "detail": detail,
                "pass": bool(ok), "required": name in NEW_KERNELS}

    rows = [
        row("G(0)=0", abs(float(respray.cdf(np.array([0.0]))[0])) <= tau,
            f"G(0)={float(respray.cdf(np.array([0.0]))[0]):.3e}"),
        row("G(1)=1", abs(float(respray.cdf(np.array([1.0]))[0]) - 1.0) <= tau,
            f"G(1)={float(respray.cdf(np.array([1.0]))[0]):.3e}"),
        row("G nondecreasing (bin edges)", bool(np.all(np.diff(G_edges) >= -tau)),
            f"min diff={float(np.diff(G_edges).min()):.3e}"),
        row("G nondecreasing (dense grid)", bool(np.all(np.diff(G_grid) >= -tau)),
            f"min diff={float(np.diff(G_grid).min()):.3e}"),
        row("q_j >= -tau_num (raw)", float(raw_q.min()) >= -tau,
            f"min q={float(raw_q.min()):.3e}"),
        row("|sum(q)-1| <= tau_num (raw)", abs(float(raw_q.sum()) - 1.0) <= tau,
            f"sum={float(raw_q.sum()):.16e}"),
        row("q_j >= -tau_num (used bins)", float(q_norm.min()) >= -tau,
            f"min q={float(q_norm.min()):.3e}"),
        row("|sum(q)-1| <= tau_num (used bins)", abs(float(q_norm.sum()) - 1.0) <= tau,
            f"sum={float(q_norm.sum()):.16e}"),
    ]
    if name == "K0":
        rows.append(row("Delta0 == 0 everywhere", float(np.abs(delta).max()) <= tau,
                        f"max|Delta|={float(np.abs(delta).max()):.3e}"))
    if name == "K-":
        rows.append(row("Delta_minus <= 0 everywhere", float(delta.max()) <= tau,
                        f"max Delta={float(delta.max()):.3e}"))
    if name == "KW":
        err = float(np.abs(delta - 0.25 * x * (1.0 - x)).max())
        rows.append(row("Delta_W(j/K) == (1/4)(j/K)(1-j/K)", err <= tau,
                        f"max err={err:.3e}"))
    if name == "KS":
        pieces = {
            "left": lambda x: 1.75 * x,
            "middle": lambda x: x + 0.25,
            "right": lambda x: 0.75 + 0.25 * x,
        }
        rows.append(row("Delta_S >= 0 everywhere", float(delta.min()) >= -tau,
                        f"min Delta={float(delta.min()):.3e}"))
        inside = (x >= 1.0 / 3.0) & (x <= 2.0 / 3.0)
        rows.append(row("flat top sampled on grid", bool(inside.any()),
                        f"grid points in [1/3, 2/3]: {int(inside.sum())}"))
        rows.append(row("max Delta_S == 1/4 on grid", bool(inside.any()) and
                        abs(float(delta[inside].max()) - 0.25) <= tau,
                        f"max Delta={float(delta[inside].max()):.16f}" if inside.any() else "n/a"))
        for bp, left_f, right_f, tag in [
            (1.0 / 3.0, pieces["left"], pieces["middle"], "1/3"),
            (2.0 / 3.0, pieces["middle"], pieces["right"], "2/3"),
        ]:
            g_left, g_right = float(left_f(bp)), float(right_f(bp))
            rows.append(row(f"G continuous at {tag}",
                            abs(g_left - g_right) <= tau,
                            f"|{g_left:.16f} - {g_right:.16f}| = {abs(g_left - g_right):.3e}"))
            rows.append(row(f"Delta_S continuous at {tag}",
                            abs((g_left - bp) - (g_right - bp)) <= tau,
                            f"|Delta_left - Delta_right| = {abs((g_left - bp) - (g_right - bp)):.3e}"))
    return rows


# --------------------- A3/A4: M1C.1-compatible integrity summary ---------------------

INTEGRITY_FIELDS = [
    "t_exact_stop", "t_U_best", "U_density_best", "U_density_at_stop",
    "improvement_retention", "moved_mass_median", "moved_mass_total",
    "D_max_p5", "D_max_p50", "D_max_p95", "U_density_p5",
    "U_density_p50", "U_density_p95", "share_delta_U_positive",
    "share_D_max_increases", "near_mass_05_at_stop",
]


def m1c1_compatible_summary(variant: str, K: int, run: dict, alpha: float) -> dict:
    """Summary fields with the exact formulas/formats of the committed
    M1C.1 variant_summary.csv (experiments/m1c1_uniform_kernel_ablation),
    so K0/KC integrity is a string comparison against committed reality."""
    stop = run["stopped_at"]
    rows = run["rows"]
    n = stop if stop is not None else len(rows) - 1
    Ud = np.array([K * r["U_L2"] for r in rows])
    Dm = np.array([r["D_max"] for r in rows[:n]])
    moved = np.array([r["removed_mass"] for r in rows[:n]])
    traj_Ud = Ud[:n + 1]
    u0, u_best = float(traj_Ud[0]), float(traj_Ud.min())
    t_best = int(traj_Ud.argmin())
    u_stop = float(traj_Ud[-1])
    Db, Ub = band(Dm), band(traj_Ud)
    last25 = slice(max(int(0.75 * n), 0), n)
    late_band = Dm < 0.005
    dU = np.diff(traj_Ud)
    return {
        "variant": variant, "K": K,
        "t_exact_stop": stop if stop is not None else "none",
        "active_rounds": n,
        "t_U_best": t_best,
        "U_density_initial": f"{u0:.6e}",
        "U_density_best": f"{u_best:.6e}",
        "U_density_at_stop": f"{u_stop:.6e}",
        "improvement_retention": f"{1.0 - (u_stop - u_best) / (u0 - u_best):.6f}",
        "D_max_p5": f"{Db['p5']:.5f}", "D_max_p50": f"{Db['p50']:.5f}",
        "D_max_p95": f"{Db['p95']:.5f}",
        "U_density_p5": f"{Ub['p5']:.5f}", "U_density_p50": f"{Ub['p50']:.5f}",
        "U_density_p95": f"{Ub['p95']:.5f}",
        "moved_mass_median": f"{float(np.median(moved)):.6f}",
        "moved_mass_median_last25pct": f"{float(np.median(moved[last25])):.6f}",
        "moved_mass_median_late_band":
            f"{float(np.median(moved[late_band])):.6f}" if late_band.any() else "n/a",
        "moved_mass_total": f"{float(moved.sum()):.6f}",
        "share_delta_U_positive": f"{float((dU > 0).mean()):.4f}",
        "share_D_max_increases": f"{float((np.diff(Dm) > 0).mean()):.4f}",
        "near_mass_05_at_stop": f"{rows[n]['near_mass_05']:.6f}",
    }


# ------------------------- per-round diagnostics + theory ----------------------------

def per_round_diagnostics(kernel: str, K: int, run: dict, g_bins: np.ndarray,
                          alpha: float, tau: float, switch_threshold: float
                          ) -> tuple[list[dict], dict]:
    """Frozen §12 diagnostics on every active round + §13 theory verification.

    R_t uses exactly the M1C.4 semantics: max over j with gamma_t(j) > 0 of
    M_t max(Delta(j), 0) / gamma_t(j), Delta from the kernel's own bins.
    The decomposition identity and the exact rebound criterion are checked
    against the ACTUAL next state; FP/FN keep plain strict semantics, with
    |worst excess| <= tau_num rounds flagged as numerical-boundary cases.
    """
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
        pos_reinj = moved * np.maximum(delta, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(gamma > 0.0,
                             pos_reinj / np.where(gamma > 0.0, gamma, 1.0), 0.0)
        R = float(ratio.max())
        j_drive = int(np.argmax(ratio)) + 1
        excess = moved * delta - gamma
        worst = float(excess.max())
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
            "worst_excess": f"{worst:.6e}",
            "numerical_boundary": bool(boundary_case),
            "j_drive": j_drive,
            "R_t": f"{R:.10e}",
            "stopped": stop is not None,
            "t_stop": stop if stop is not None else "",
        })
    return out, {
        "active_rounds": n, "TP": tp, "TN": tn, "FP": fp, "FN": fn,
        "numerical_boundary": n_boundary, "FP_nonboundary": fp_nb,
        "FN_nonboundary": fn_nb, "max_decomposition_residual": max_resid,
        "max_mass_error": max_mass_err, "min_bin_mass": min_bin,
    }


def classify(N_R: int, S_R_over_K: float | None, cfg: dict) -> str:
    """Frozen descriptive classification (config.classification_convention)."""
    if N_R == 0:
        return "No rebound"
    if N_R == 1:
        return "Isolated rebound"
    if S_R_over_K is not None and S_R_over_K >= 1.0:
        return "Long-span recurrent / churn-like transient"
    return "Short recurrent transient"


def trajectory_summary(kernel: str, K: int, run: dict, per_rows: list[dict],
                       theory: dict, cfg: dict) -> dict:
    """Rebound / recurrence summary (§15) + Table 2 row (§21)."""
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
    n_switch = sum(1 for s in switches if s)
    moved = np.array([float(r["moved_mass"]) for r in per_rows])
    d0 = np.array([float(r["D_max_before"]) for r in per_rows])
    # longest contraction-only run (consecutive non-rebound active rounds)
    longest_run = cur = 0
    for r in per_rows:
        cur = 0 if r["rebound"] else cur + 1
        longest_run = max(longest_run, cur)
    horizon = 10 * K
    return {
        "kernel": kernel, "K": K,
        "active_rounds": n,
        "stop_status": ("exact_stop" if stop is not None
                        else "not stopped within pre-registered Gate-1 horizon"),
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
        "boundary_switch_count": n_switch,
        "boundary_switch_freq": f"{n_switch / max(len(switches), 1):.5f}",
        "abs_delta_a_median": f"{float(np.median(da)):.5f}" if len(da) else "",
        "abs_delta_a_p95": f"{float(np.percentile(da, 95)):.5f}" if len(da) else "",
        "moved_mass_median": f"{float(np.median(moved)):.6f}",
        "moved_mass_total": f"{float(moved.sum()):.6f}",
        "D_max_initial": f"{d0[0]:.6f}" if len(d0) else "",
        "D_max_max": f"{d0.max():.6f}" if len(d0) else "",
        "D_max_final": f"{run['rows'][n]['D_max']:.6e}",
        "max_decomposition_residual": f"{theory['max_decomposition_residual']:.3e}",
        "criterion_TP": theory["TP"], "criterion_TN": theory["TN"],
        "criterion_FP": theory["FP"], "criterion_FN": theory["FN"],
        "numerical_boundary_count": theory["numerical_boundary"],
        "FP_nonboundary": theory["FP_nonboundary"],
        "FN_nonboundary": theory["FN_nonboundary"],
        "max_mass_error": f"{theory['max_mass_error']:.3e}",
        "min_bin_mass": f"{theory['min_bin_mass']:.3e}",
        "qualitative_category": classify(n_reb, S_R_over_K, cfg),
    }


# ------------------------------------ figures ----------------------------------------

def _panels(kernels: list[str], K: int, results: dict):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2))
    flat = np.asarray(axes).reshape(-1)
    for ax, name in zip(flat, kernels):
        res = results[(name, K)]
        yield ax, name, res


def _traj_arrays(res: dict) -> dict:
    run, per = res["run"], res["per"]
    n = res["theory"]["active_rounds"]
    rows = run["rows"]
    d_traj = np.array([rows[t]["D_max"] for t in range(n + 1)])
    a = np.array([float(r["a_t"]) for r in per])
    R = np.array([float(r["R_t"]) for r in per])
    reb = np.array([bool(r["rebound"]) for r in per])
    return {"d": d_traj, "a": a, "R": R, "reb": reb, "n": n}


def plot_fig_D_max(results: dict, K: int, kernels: list[str], path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), sharex=False)
    for ax, name in zip(np.asarray(axes).reshape(-1), kernels):
        arr = _traj_arrays(results[(name, K)])
        ax.semilogy(np.arange(arr["n"] + 1), np.maximum(arr["d"], 1e-16),
                    color=COLORS[name], lw=0.9)
        ax.set_title(f"{name}  (K={K})", fontsize=10)
        ax.set_ylabel("D_max")
        ax.set_xlabel("round t")
    fig.suptitle(f"M2 Gate 1 Figure A — D_max(t) at K={K} "
                 "(semilogy; K0 omitted: monotone drop to ~1e-13 in 6 rounds)")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_fig_boundary(results: dict, K: int, kernels: list[str], path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), sharex=False)
    for ax, name in zip(np.asarray(axes).reshape(-1), kernels):
        arr = _traj_arrays(results[(name, K)])
        ax.plot(np.arange(arr["n"]), arr["a"], color=COLORS[name],
                lw=0.5, marker=".", ms=0.6, ls="None" if arr["n"] > 400 else "-")
        ax.set_title(f"{name}  (K={K})", fontsize=10)
        ax.set_ylabel("boundary a_t")
        ax.set_xlabel("round t")
        ax.set_ylim(0, 1)
    fig.suptitle(f"M2 Gate 1 Figure B — sweep boundary a_t at K={K}")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_fig_rebound(results: dict, K: int, kernels: list[str], path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), sharex=False)
    for ax, name in zip(np.asarray(axes).reshape(-1), kernels):
        arr = _traj_arrays(results[(name, K)])
        t = np.arange(arr["n"])
        ax.semilogy(t, np.maximum(arr["R"], 1e-3), color=COLORS[name],
                    lw=0.5, marker=".", ms=0.6, ls="None" if arr["n"] > 400 else "-")
        reb_mask = arr["reb"]
        if reb_mask.any():
            ax.scatter(t[reb_mask], np.maximum(arr["R"][reb_mask], 1e-3),
                       s=4, color="black", zorder=3, label="rebound rounds")
            ax.legend(fontsize=7, loc="best")
        ax.axhline(1.0, color="gray", lw=0.8, ls="--")
        f_R = float(arr["reb"].mean()) if arr["n"] else 0.0
        ax.set_title(f"{name}  (K={K}, N_R={int(arr['reb'].sum())}, "
                     f"f_R={f_R:.3f})", fontsize=10)
        ax.set_ylabel("R_t")
        ax.set_xlabel("round t")
    fig.suptitle(f"M2 Gate 1 Figure C — rebound pressure R_t and rebound markers at K={K}")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_fig_comparison(results: dict, K: int, kernels: list[str], path: Path) -> None:
    """KC vs KW vs KS (plus K- as reference) direct comparison at K=100."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))
    for name in kernels:
        arr = _traj_arrays(results[(name, K)])
        ax1.semilogy(np.arange(arr["n"] + 1), np.maximum(arr["d"], 1e-16),
                     color=COLORS[name], lw=0.9, label=name)
        ax2.plot(np.arange(arr["n"]), arr["a"], color=COLORS[name],
                 lw=0.4, label=name)
    ax1.set_xlabel("round t")
    ax1.set_ylabel("D_max")
    ax1.set_yscale("log")
    ax1.legend(fontsize=8)
    ax2.set_xlabel("round t")
    ax2.set_ylabel("boundary a_t")
    ax2.set_ylim(0, 1)
    ax2.legend(fontsize=8)
    fig.suptitle(f"M2 Gate 1 — KC vs KW vs KS (K- reference) at K={K}: "
                 "D_max (left) and boundary (right)")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------- main ---------------------------------------

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
    registry = build_registry(cfg)
    started = time.time()

    print(f"M2 Gate 1 redistribution-mismatch generality — config: {config_path}")
    print(f"  alpha={alpha}, exact tol={tol}, tau_num={tau}, "
          f"horizon H(K)=10K, primary K={cfg['primary_K']}")

    results: dict[tuple[str, int], dict] = {}
    per_round_all: list[dict] = []
    summary_rows: list[dict] = []

    def run_and_diagnose(kernel: str, K: int) -> dict:
        respray = registry[kernel]["respray"]
        g_bins = (respray if respray is not None else q_near).bin_probs(K)
        run = run_m1a_deterministic(
            q_near, K=K, T_max=10 * K, alpha=alpha, stop_tolerance=tol,
            redistribution="throw",
            redistribution_dist=respray,
        )
        per, theory = per_round_diagnostics(kernel, K, run, g_bins, alpha,
                                            tau, switch_threshold)
        summary = trajectory_summary(kernel, K, run, per, theory, cfg)
        per_round_all.extend(per)
        summary_rows.append(summary)
        stop = run["stopped_at"]
        print(f"  [{kernel}] K={K}: active={theory['active_rounds']}, "
              f"stop={'t=' + str(stop) if stop is not None else 'CENSORED at 10K'}, "
              f"N_R={summary['N_R']}, f_R={summary['f_R']}, "
              f"resid={summary['max_decomposition_residual']}, "
              f"FP={theory['FP']}, FN={theory['FN']} "
              f"(non-boundary {theory['FP_nonboundary']}/{theory['FN_nonboundary']}), "
              f"category={summary['qualitative_category']}")
        entry = {"run": run, "per": per, "theory": theory, "summary": summary}
        results[(kernel, K)] = entry
        return entry

    # ---------------- A2: static kernel validation ----------------
    print("\n=== A2: static kernel validation (tau_num = %.0e) ===" % tau)
    static_rows: list[dict] = []
    for name in NEW_KERNELS + ["K0", "KC"]:
        respray = registry[name]["respray"]
        respray = respray if respray is not None else q_near
        for K in registry[name]["K_list"]:
            static_rows.extend(static_kernel_checks(name, respray, K, tau))
    write_csv(out_dir / "static_kernel_validation.csv",
              ["kernel", "K", "check", "detail", "pass", "required"], static_rows)
    required_fail = [r for r in static_rows if r["required"] and not r["pass"]]
    n_req = sum(1 for r in static_rows if r["required"])
    print(f"  static checks: {n_req - len(required_fail)}/{n_req} required passed "
          f"({len(static_rows) - len(required_fail)}/{len(static_rows)} incl. reference rows)")
    if required_fail:
        for r in required_fail:
            print(f"  [FAIL] {r['kernel']} K={r['K']}: {r['check']} ({r['detail']})")
        print("STATIC VALIDATION FAILED — no dynamics will be run.")
        return 1

    # ---------------- A3/A4: K0/KC baseline integrity ----------------
    print("\n=== A3/A4: K0/KC baseline integrity vs committed M1C.1/M1B.2 ===")
    for kernel in ["K0", "KC"]:
        for K in CANONICAL_K_LIST:
            run_and_diagnose(kernel, K)
    committed = {(r["variant"], int(r["K"])): r for r in read_csv_rows(
        REPO_ROOT / cfg["baseline_integrity_sources"]["m1c1_variant_summary"])}
    ann = {int(r["K"]): r for r in read_csv_rows(
        REPO_ROOT / cfg["baseline_integrity_sources"]["m1b2_annihilation_events"])}
    integ_rows: list[dict] = []
    variant_of = {"K0": "U", "KC": "A"}
    for kernel, variant in variant_of.items():
        for K in CANONICAL_K_LIST:
            fresh = m1c1_compatible_summary(variant, K, results[(kernel, K)]["run"],
                                            alpha)
            for f in INTEGRITY_FIELDS:
                exp = str(committed[(variant, K)][f])
                got = str(fresh[f])
                integ_rows.append({"check": f"{kernel} ({variant}) K={K}: {f}",
                                   "expected": exp, "got": got,
                                   "pass": exp == got})
    for K in CANONICAL_K_LIST:
        stop_kc = results[("KC", K)]["run"]["stopped_at"]
        integ_rows.append({"check": f"KC K={K}: t_exact_stop vs M1B.2",
                           "expected": ann[K]["t_stop"], "got": str(stop_kc),
                           "pass": str(stop_kc) == ann[K]["t_stop"]})
    n_pass = sum(1 for r in integ_rows if r["pass"])
    print(f"  baseline integrity: {n_pass}/{len(integ_rows)} checks passed")
    write_csv(out_dir / "baseline_integrity.csv",
              ["check", "expected", "got", "pass"], integ_rows)
    if n_pass < len(integ_rows):
        for r in integ_rows:
            if not r["pass"]:
                print(f"  [FAIL] {r['check']}: expected {r['expected']}, got {r['got']}")
        print("BASELINE DRIFT — Gate 1 scientific execution stops here "
              "(config.baseline_policy). Locate the regression first.")
        return 1

    # ---------------- A5: K- negative control at K=100 ----------------
    print("\n=== A5: K- negative control @ K=100 (theorem-backed acceptance gate) ===")
    km = run_and_diagnose("K-", 100)
    th, summ = km["theory"], km["summary"]
    strict_ok = all(float(r["d_D_max"]) < 0.0 for r in km["per"])
    min_decrease = min(float(r["d_D_max"]) for r in km["per"])
    checks = {
        "zero genuine rebound": summ["N_R"] == 0,
        "strict active-round D_max decrease": strict_ok,
        "decomposition residual <= tau_num":
            th["max_decomposition_residual"] <= tau,
        "rebound criterion FP=FN=0": th["FP"] == 0 and th["FN"] == 0,
        "mass conservation <= tau_num": th["max_mass_error"] <= tau,
    }
    for label, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    print(f"  (min per-round D_max decrease: {min_decrease:.3e})")
    if not all(checks.values()):
        print("K- ACCEPTANCE FAILED — do NOT interpret KW/KS scientific results; "
              "debug implementation / theory helpers first (config.acceptance_rule_K_minus).")
        write_csv(out_dir / "per_round_diagnostics.csv",
                  list(per_round_all[0]), per_round_all)
        write_csv(out_dir / "rebound_summary.csv",
                  list(summary_rows[0]), summary_rows)
        return 1

    # ---------------- B1/B2 + C: KW / KS ----------------
    print("\n=== B1/B2: KW, KS @ K=100 (main screening) ===")
    for kernel in ["KW", "KS"]:
        run_and_diagnose(kernel, 100)
    print("\n=== C: KW, KS cross-K confirmation @ K = 50, 200, 400 ===")
    for K in [50, 200, 400]:
        for kernel in ["KW", "KS"]:
            run_and_diagnose(kernel, K)

    _write_outputs(out_dir, cfg, per_round_all, summary_rows, results, started)

    print("\n=== Table 2 (dynamics summary) ===")
    for r in summary_rows:
        print(f"  {r['kernel']:>3} K={r['K']:>3}: {r['stop_status'][:24]:<24} "
              f"t_stop={str(r['t_stop']):>5} N_R={r['N_R']:>4} f_R={r['f_R']} "
              f"S_R/K={str(r['S_R_over_K']):>8} sw={r['boundary_switch_freq']} "
              f"cat={r['qualitative_category']}")

    print(f"\noutputs written to: {out_dir}")
    print(f"runtime: {round(time.time() - started, 2)} s")
    return 0


def _write_outputs(out_dir: Path, cfg: dict, per_round_all: list[dict],
                   summary_rows: list[dict], results: dict, started: float) -> None:
    write_csv(out_dir / "per_round_diagnostics.csv",
              list(per_round_all[0]), per_round_all)
    write_csv(out_dir / "rebound_summary.csv",
              list(summary_rows[0]), summary_rows)
    K_primary = int(cfg["primary_K"])
    plot_fig_D_max(results, K_primary, ["KC", "K-", "KW", "KS"],
                   out_dir / "fig_gate1_D_max_K100.png")
    plot_fig_boundary(results, K_primary, ["KC", "K-", "KW", "KS"],
                      out_dir / "fig_gate1_boundary_K100.png")
    plot_fig_rebound(results, K_primary, ["KC", "K-", "KW", "KS"],
                     out_dir / "fig_gate1_rebound_R_K100.png")
    plot_fig_comparison(results, K_primary, ["KC", "K-", "KW", "KS"],
                        out_dir / "fig_gate1_KC_KW_KS_comparison_K100.png")

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
            "dynamics": {f"{r['kernel']}@K={r['K']}": {
                "stop_status": r["stop_status"], "t_stop": r["t_stop"],
                "N_R": r["N_R"], "f_R": r["f_R"], "S_R_over_K": r["S_R_over_K"],
                "category": r["qualitative_category"]}
                for r in summary_rows},
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "git": git_state(),
        "runtime_seconds": round(time.time() - started, 2),
        "command": "python experiments/m2_gate1_generality/run_experiment.py",
        "outputs": {name: sha256_of(out_dir / name) for name in output_files},
    }
    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )


if __name__ == "__main__":
    raise SystemExit(main())
