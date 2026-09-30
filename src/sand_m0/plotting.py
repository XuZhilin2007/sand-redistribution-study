"""M0 figures: per-experiment uniformity curves and deterministic-vs-MC panels."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless rendering only; no interactive backend needed

import matplotlib.pyplot as plt
import numpy as np

from .mechanism import quadratic_coefficients

# Values at or below the floor (e.g. deterministic U = 0 for the uniform
# control at t = 0) are drawn at the axis floor on log-scale plots.
_FLOOR = 1e-10


def _safe_log(values: np.ndarray) -> np.ndarray:
    return np.maximum(np.asarray(values, dtype=float), _FLOOR)


def plot_experiment(
    det_rows: list[dict],
    mc_mean: np.ndarray,
    mc_std: np.ndarray,
    *,
    exp_id: str,
    dist_label: str,
    noise_floor: float,
    det_best: int,
    mc_best: int,
    path: Path,
) -> None:
    """Uniformity error vs round: deterministic line, MC mean and +-1 sd band."""
    ts = np.array([row["t"] for row in det_rows])
    det_u = _safe_log([row["U_L2"] for row in det_rows])

    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=150)
    lower = _safe_log(mc_mean - mc_std)
    upper = _safe_log(mc_mean + mc_std)
    ax.fill_between(ts, lower, upper, alpha=0.25, linewidth=0, label="MC ±1 sd (across seeds)")
    ax.plot(ts, _safe_log(mc_mean), "o-", ms=3.5, lw=1.4, color="tab:blue", label="Monte Carlo mean")
    ax.plot(ts, det_u, "-", lw=1.8, color="black", label="Deterministic (expected mass)")
    ax.axhline(noise_floor, color="gray", ls=":", lw=1.0,
               label=f"finite-N noise floor ≈ {noise_floor:.2e}")

    ax.axvline(det_best, color="black", ls="--", lw=0.8)
    ax.axvline(mc_best, color="tab:blue", ls="--", lw=0.8)
    ymin, ymax = ax.get_ylim()
    ax.annotate(f"det best t={det_best}", xy=(det_best, ymin), xytext=(3, 12),
                textcoords="offset points", fontsize=8, rotation=90, va="bottom")
    ax.annotate(f"MC mean best t={mc_best}", xy=(mc_best, ymin), xytext=(3, 12),
                textcoords="offset points", fontsize=8, rotation=90, va="bottom",
                color="tab:blue")

    ax.set_yscale("log")
    ax.set_xlabel("round t")
    ax.set_ylabel(r"$U(t)=\sum_i (p_i - 1/K)^2$   (log scale)")
    ax.set_title(f"{exp_id} — {dist_label}\nuniformity error vs round (M0 baseline)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_det_vs_mc(
    curves: dict[str, dict],
    *,
    path: Path,
    N: int,
    n_seeds: int,
) -> None:
    """One panel per experiment: deterministic trajectory vs MC mean.

    `curves` maps exp_id -> {"label", "det_u", "mc_mean"}.
    """
    exp_ids = sorted(curves)
    fig, axes = plt.subplots(1, len(exp_ids), figsize=(4.6 * len(exp_ids), 4.4),
                             dpi=150, sharey=False)
    if len(exp_ids) == 1:
        axes = [axes]
    for ax, exp_id in zip(axes, exp_ids):
        c = curves[exp_id]
        ax.plot(c["det_u"], "-", lw=1.8, color="black", label="Deterministic")
        ax.plot(c["mc_mean"], "o-", ms=3, lw=1.2, color="tab:blue", label="MC mean")
        ax.set_title(f"{exp_id}\n{c['label']}", fontsize=9)
        ax.set_xlabel("round t")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    axes[0].set_ylabel(r"$U(t)=\sum_i (p_i - 1/K)^2$")
    fig.suptitle(f"Deterministic vs Monte Carlo mean (N={N}, {n_seeds} seeds) — "
                 "curves nearly coincide, differences are finite-sample noise",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path)
    plt.close(fig)


# ------------------------- M0.1 mechanism figures -------------------------

_ALPHA_COLORS = {
    0.0: "lightgray",
    0.10: "tab:green",
    0.25: "tab:blue",
    0.50: "tab:orange",
    0.75: "tab:red",
    1.00: "tab:purple",
}


def plot_alpha_rounds(
    det_curves: dict[float, dict], *, path: Path, T: int
) -> None:
    """U vs round t for every alpha (M0.1 time-scale view, log y)."""
    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=150)
    for alpha, c in sorted(det_curves.items()):
        color = _ALPHA_COLORS.get(alpha, None)
        ls = "--" if alpha == 0.0 else "-"
        lw = 1.0 if alpha == 0.0 else 1.6
        ax.plot(range(len(c["u"])), c["u"], ls, lw=lw, color=color,
                marker="o" if alpha != 0.0 else None, ms=3,
                label=fr"$\alpha={alpha:g}$")
        if alpha != 0.0:
            ax.axvline(c["best_t"], color=color, ls=":", lw=0.8)
            ax.annotate(f"t={c['best_t']}", xy=(c["best_t"], ax.get_ylim()[0]),
                        xytext=(2, 10), textcoords="offset points",
                        fontsize=7, rotation=90, color=color, va="bottom")
    ax.set_yscale("log")
    ax.set_xlabel("round t")
    ax.set_ylabel(r"$U(t)=\sum_i (p_i - 1/K)^2$   (log scale)")
    ax.set_title("M0.1 — correction strength changes the pace, not the shape\n"
                 r"deterministic $U(t)$ for $\alpha\in\{0,0.1,0.25,0.5,0.75,1\}$"
                 "  (q_near, a=0.5, K=100)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_rescaling(
    det_curves: dict[float, dict],
    *,
    z_grid: np.ndarray,
    u_disc: list[float],
    u_cont: list[float],
    z_star: float,
    constants,
    path: Path,
) -> None:
    """U vs state variable z: all alphas collapse onto one exact parabola."""
    fig, (ax, axr) = plt.subplots(1, 2, figsize=(12.0, 5.0), dpi=150,
                                  gridspec_kw={"width_ratios": (3, 2)})
    ax.plot(z_grid, u_disc, "-", lw=2.0, color="black",
            label="exact discrete quadratic $U(z)$")
    ax.plot(z_grid, u_cont, "--", lw=1.2, color="gray",
            label="continuous limit ($z^*=3/4$ exactly)")
    for alpha, c in sorted(det_curves.items()):
        if alpha == 0.0:
            continue
        ax.plot(c["z"], c["u"], "o", ms=4, mfc="none",
                color=_ALPHA_COLORS.get(alpha), alpha=0.85,
                label=fr"simulated $\alpha={alpha:g}$")
    ax.axvline(0.75, color="tab:red", ls=":", lw=1.0, label=r"$z^*=3/4$ (continuous)")
    ax.axvline(z_star, color="tab:red", ls="-.", lw=1.0,
               label=f"discrete $z^*={z_star:.6f}$")
    ax.axvline(0.75 ** 4, color="tab:purple", ls=":", lw=0.8)
    ax.axvline(0.75 ** 5, color="tab:purple", ls=":", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel(r"state variable $z=\rho^t$ (near-zone mass $=0.75\,z$)")
    ax.set_ylabel(r"$U(z)=\sum_i (p_i - 1/K)^2$   (log scale)")
    ax.set_title("M0.1 — time rescaling\nall $\\alpha$ on one exact state path", fontsize=10)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7, loc="upper right")

    A, B, Cc = quadratic_coefficients(constants)
    for alpha, c in sorted(det_curves.items()):
        if alpha == 0.0:
            continue
        resid = [abs(u - (A * z * z + B * z + Cc)) for u, z in zip(c["u"], c["z"])]
        axr.plot(c["z"], resid, "o", ms=3.5, mfc="none", alpha=0.85,
                 color=_ALPHA_COLORS.get(alpha), label=fr"$\alpha={alpha:g}$")
    axr.set_yscale("log")
    axr.set_xlabel(r"state variable $z$")
    axr.set_ylabel(r"$|U_{\rm sim} - U(z)|$")
    axr.set_title("residual vs exact quadratic\n(floating-point level)", fontsize=10)
    axr.grid(True, which="both", alpha=0.25)
    axr.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_mc_spotcheck(
    *,
    mc_csv: Path,
    z_grid: np.ndarray,
    u_disc: list[float],
    noise_floor: float,
    path: Path,
) -> None:
    """MC spot checks against the exact state-path curve (log y)."""
    import csv as _csv

    with mc_csv.open(encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh))
    by_alpha: dict[float, list[tuple[float, float]]] = {}
    for r in rows:
        by_alpha.setdefault(float(r["alpha"]), []).append(
            (float(r["z"]), float(r["U_L2"]))
        )

    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=150)
    ax.plot(z_grid, u_disc, "-", lw=2.0, color="black",
            label="exact discrete quadratic $U(z)$")
    for alpha, pts in sorted(by_alpha.items()):
        zs, us = zip(*pts)
        ax.plot(zs, us, "o", ms=3, mfc="none", alpha=0.7,
                color=_ALPHA_COLORS.get(alpha), label=fr"MC $\alpha={alpha:g}$ (5 seeds)")
    ax.axhline(noise_floor, color="gray", ls=":", lw=1.0,
               label=f"finite-N noise floor ≈ {noise_floor:.2e}")
    ax.set_yscale("log")
    ax.set_xlabel(r"state variable $z=\rho^t$")
    ax.set_ylabel(r"$U=\sum_i (p_i - 1/K)^2$   (log scale)")
    ax.set_title("M0.1 — Monte Carlo spot checks on the state path\n"
                 r"(q_near, N=100000, seeds 20260917–21)", fontsize=10)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------- M0.2 Phase 2A family figures ---------------------

_FAMILY_COLORS = {
    -0.25: "tab:green",
    0.0: "gray",
    0.25: "tab:cyan",
    0.5: "tab:orange",
    1.0: "tab:blue",
    2.0: "tab:red",
}


def _family_color(p: float) -> str:
    return _FAMILY_COLORS.get(p, None)


def plot_family_shapes(dists: dict[float, object], *, a: float, path: Path) -> None:
    """Initial density shapes of every frozen family member."""
    xs = np.linspace(0.0, 1.0, 1001)
    fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=150)
    ax.axvspan(0.0, a, color="tab:blue", alpha=0.08,
               label=f"near zone x < {a}")
    for p, dist in sorted(dists.items()):
        cdf = dist.cdf
        dens = np.diff(cdf(xs)) / (xs[1] - xs[0])
        ax.plot(xs[:-1], dens, lw=1.8, color=_family_color(p),
                label=fr"$p={p:g}$")
    ax.set_xlabel("x  (0 = nearest to operator)")
    ax.set_ylabel(r"density of $q_p(x) = (p+1)(1-x)^p$")
    ax.set_title("M0.2 — frozen power family: initial throw distributions")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_u_vs_round_family(
    curves: dict[float, dict], *, alpha: float, path: Path
) -> None:
    """U vs round for every family member (baseline alpha), log y."""
    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=150)
    for p, c in sorted(curves.items()):
        ax.plot(range(len(c["u"])), c["u"], "-o", ms=3, lw=1.5,
                color=_family_color(p), label=fr"$p={p:g}$ (best t={c['best_t']})")
        ax.axvline(c["best_t"], color=_family_color(p), ls=":", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("round t")
    ax.set_ylabel(r"$U(t)=\sum_i (p_i - 1/K)^2$   (log scale)")
    ax.set_title(fr"M0.2 — uniformity vs round, power family ($\alpha={alpha:g}$, "
                 "a=0.5, K=100)\ninterior optimum only when near-zone shape is more "
                 "concentrated per unit mass (p > 0)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_u_vs_z_family(
    curves: dict[float, dict],
    *,
    theory: dict[float, tuple[float, float, float]],
    alpha_trajectories: dict[float, dict[float, tuple[list[float], list[float]]]],
    path: Path,
) -> None:
    """U vs state variable z: per-q exact quadratics; alpha points collapse.

    `curves` maps p -> {"z", "u"} for the baseline alpha. `theory` maps
    p -> (A, B, C) of that q's exact quadratic U(z). `alpha_trajectories`
    maps p -> {alpha: (zs, us)} for the spot-checked alphas.
    """
    fig, ax = plt.subplots(figsize=(9.5, 5.6), dpi=150)
    for p, (A, B, C) in sorted(theory.items()):
        zs = np.linspace(0.0, 1.0, 400)
        ax.plot(zs, [A * z * z + B * z + C for z in zs], "-", lw=1.3,
                color=_family_color(p), alpha=0.8, label=fr"$U_{{p={p:g}}}(z)$ theory")
        c = curves[p]
        ax.plot(c["z"], c["u"], "o", ms=4, mfc="none", color=_family_color(p),
                label=fr"sim $p={p:g}$ ($\alpha=0.25$)")
        for alpha, (zs_a, us_a) in sorted(alpha_trajectories.get(p, {}).items()):
            ax.plot(zs_a, us_a, "s", ms=3.5, mfc="none", alpha=0.65,
                    color=_family_color(p),
                    label=fr"sim $p={p:g}$ ($\alpha={alpha:g}$)")
    ax.set_yscale("log")
    ax.set_xlabel(r"state variable $z=\rho^t$")
    ax.set_ylabel(r"$U(z)$   (log scale)")
    ax.set_title("M0.2 — one exact quadratic per distribution; every alpha moves\n"
                 "points along the same curve (time rescaling, spot-checked)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=6.5, ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_prediction_vs_observed(
    comp_rows: list[dict], alpha_rows: list[dict], *, path: Path
) -> None:
    """Identity panels: predicted vs observed best round (and z* where interior)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.8), dpi=150)
    lims = [-1, 13]
    ax1.plot(lims, lims, "k--", lw=1.0, label="perfect prediction")
    for row in comp_rows:
        p = float(row["p"])
        ax1.plot(row["predicted_best_round"], row["observed_best_round"], "o",
                 ms=9, mfc="none", color=_family_color(p))
        ax1.annotate(fr"p={p:g}", xy=(row["predicted_best_round"],
                     row["observed_best_round"]), xytext=(5, 5),
                     textcoords="offset points", fontsize=8, color=_family_color(p))
    for row in alpha_rows:
        p = float(row["p"])
        if row["t_star_continuous"] == "":
            continue
        ax1.plot(round(float(row["t_star_continuous"])), row["observed_best_round"],
                 "s", ms=7, mfc="none", color=_family_color(p), alpha=0.6)
    ax1.set_xlim(lims)
    ax1.set_ylim(lims)
    ax1.set_xlabel("predicted best round (theory, before simulation)")
    ax1.set_ylabel("observed best round (canonical simulator)")
    ax1.set_title("predicted vs observed best round\n"
                  "(circles: baseline alpha; squares: alpha spot checks)", fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    for row in comp_rows:
        p = float(row["p"])
        pz = float(row["predicted_z_star"])
        oz = float(row["observed_best_state_z"])
        ax2.plot(pz, oz, "o", ms=9, mfc="none", color=_family_color(p))
        ax2.annotate(fr"p={p:g}", xy=(pz, oz), xytext=(5, 5),
                     textcoords="offset points", fontsize=8, color=_family_color(p))
    lims2 = [0.6, 1.4]
    ax2.plot(lims2, lims2, "k--", lw=1.0, label="perfect prediction")
    ax2.axvline(1.0, color="gray", ls=":", lw=1.0)
    ax2.axhline(1.0, color="gray", ls=":", lw=1.0)
    ax2.set_xlim(lims2)
    ax2.set_ylim(lims2)
    ax2.set_xlabel(r"predicted $z^*$ (theory)")
    ax2.set_ylabel(r"observed best-state $z$")
    ax2.set_title(r"predicted vs observed optimal state position"
                  "\n" r"(z>1: no interior optimum, t=0 best)", fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------- M1A adaptive-boundary figures ---------------------


def plot_m0_vs_m1a(
    m0_u: list[float],
    m1a_u: list[float],
    m0_best: int,
    m1a_best: int,
    *,
    path: Path,
) -> None:
    """Uniformity vs round: fixed-boundary M0 against adaptive M1A (log y)."""
    fig, ax = plt.subplots(figsize=(8.5, 5.2), dpi=150)
    ax.plot(range(len(m0_u)), m0_u, "-", lw=1.7, color="tab:red",
            label=f"M0 fixed a=0.5 (best t={m0_best}, then depletes)")
    ax.plot(range(len(m1a_u)), m1a_u, "-", lw=1.4, color="tab:blue",
            label=f"M1A adaptive boundary (best t={m1a_best})")
    ax.axvline(m0_best, color="tab:red", ls=":", lw=0.9)
    ax.axvline(m1a_best, color="tab:blue", ls=":", lw=0.9)
    ax.set_yscale("log")
    ax.set_xlabel("round t")
    ax.set_ylabel(r"$U(t)=\sum_i (p_i - 1/K)^2$   (log scale)")
    ax.set_title("M1A vs M0 — adaptive boundary removes the overshoot\n"
                 r"(q_near, $\alpha=0.25$, K=100, T=100; same metric and code)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_boundary_and_excess(
    a_seqs: dict[float, list[float]],
    u_seqs: dict[float, list[float]],
    d_max: dict[int, float],
    *,
    alpha0: float,
    tol: float,
    path: Path,
) -> None:
    """Boundary trajectory a_t (per alpha) and baseline D_max(t)."""
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8.5, 7.6), dpi=150, sharex=True)
    for alpha, seq in sorted(a_seqs.items()):
        ax.plot(range(len(seq)), seq, lw=0.9, alpha=0.85,
                label=fr"$\alpha={alpha:g}$")
    ax.set_ylabel(r"sweep boundary $a_t$")
    ax.set_ylim(0, 1.02)
    ax.set_title("M1A — the boundary follows the current state; alpha changes the "
                 "trajectory itself\n"
                 r"(first decision identical at $a=0.5$; sequences diverge from the "
                 "second decision onward)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="upper right")

    ts = sorted(d_max)
    ax2.plot(ts, [d_max[t] for t in ts], lw=1.4, color="tab:blue",
             label=fr"baseline $\alpha={alpha0:g}$")
    ax2.axhline(tol, color="gray", ls=":", lw=1.0,
                label=f"stop tolerance {tol:.0e}")
    ax2.set_yscale("log")
    ax2.set_xlabel("round t")
    ax2.set_ylabel(r"$D_{\max}(t)=\max_j D_t(j)$   (log scale)")
    ax2.set_title("cumulative near-side excess stays positive: no natural stop "
                  "within T_max", fontsize=10)
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_initial_excess(
    D_curves: dict[str, np.ndarray],
    *,
    K: int,
    path: Path,
) -> None:
    """Initial cumulative excess D(j) for the three distributions."""
    xs = np.arange(1, K + 1) / K
    fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=150)
    colors = {"near": "tab:blue", "uniform": "gray", "far": "tab:green"}
    for name, D in D_curves.items():
        ax.plot(xs, D, lw=1.8, color=colors.get(name),
                label=f"q_{name}" + ("  (stops immediately)" if name != "near"
                                     else "  (excess peaks at a=0.5)"))
    ax.axhline(0.0, color="black", lw=0.9)
    ax.set_xlabel("prefix end x = j/K")
    ax.set_ylabel(r"$D(j)=$ cumulative mass $- j/K$")
    ax.set_title("M1A — initial cumulative excess: why the controls stop at t=0\n"
                 "sweep active only while some prefix holds more than its uniform share")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------- M1A.1 long-horizon diagnostic figures ---------------------


def plot_long_horizon(
    U: np.ndarray,
    D_max: np.ndarray,
    a_t: np.ndarray,
    *,
    tol: float,
    path: Path,
    checkpoints: list[int] | None = None,
) -> None:
    """Long-horizon M1A: U (log), D_max (log) and the boundary trajectory."""
    t = np.arange(len(U))
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.6), dpi=150, sharex=True)
    axes[0].plot(t, U, lw=0.8, color="tab:blue")
    axes[0].set_yscale("log")
    axes[0].set_ylabel(r"$U(t)$ (log)")
    axes[0].set_title("M1A.1 — long-horizon deterministic diagnostics (q_near, "
                      r"$\alpha=0.25$, K=100, T=5000)")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[1].plot(t, D_max, lw=0.8, color="tab:orange")
    axes[1].axhline(tol, color="gray", ls=":", lw=1.0, label=f"stop tolerance {tol:.0e}")
    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"$D_{\max}(t)$ (log)")
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend(fontsize=8)
    axes[2].plot(t, a_t, lw=0.5, color="tab:green")
    axes[2].set_ylim(0, 1.02)
    axes[2].set_ylabel(r"boundary $a_t$")
    axes[2].set_xlabel("round t")
    axes[2].grid(True, alpha=0.3)
    if checkpoints:
        for ax in axes:
            for cp in checkpoints:
                ax.axvline(cp, color="gray", ls=":", lw=0.5, alpha=0.6)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_k_sensitivity(runs: dict[int, dict], K_list: list[int], *, path: Path) -> None:
    """Cross-resolution comparison: U_density and D_max trajectories."""
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(9.0, 7.4), dpi=150, sharex=True)
    colors = {50: "tab:green", 100: "tab:blue", 200: "tab:orange", 400: "tab:red"}
    for K in sorted(K_list):
        rows = runs[K]["rows"]
        t = np.arange(len(rows))
        ax.plot(t, K * np.array([r["U_L2"] for r in rows]), lw=0.8,
                color=colors.get(K), label=f"K={K}")
        ax2.plot(t, [r["D_max"] for r in rows], lw=0.8,
                 color=colors.get(K), label=f"K={K}")
    ax.set_yscale("log")
    ax.set_ylabel(r"$U_{\rm density}(t)=K\,U(t)$ (log)")
    ax.set_title("M1A.1 — resolution sensitivity (T=5000, same continuous q and rule)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    ax2.set_yscale("log")
    ax2.set_ylabel(r"$D_{\max}(t)$ (log)")
    ax2.set_xlabel("round t")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_d_profiles(D_curves: dict[str, np.ndarray], *, K: int, path: Path) -> None:
    """A few representative cumulative-excess profiles D(x)."""
    xs = np.arange(1, K + 1) / K
    n = len(D_curves)
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 3.4), dpi=150, sharey=False)
    if n == 1:
        axes = [axes]
    for ax, (label, D) in zip(axes, D_curves.items()):
        ax.plot(xs, D, lw=1.3, color="tab:blue")
        j = int(np.argmax(D))
        ax.axvline((j + 1) / K, color="tab:red", ls=":", lw=1.0)
        ax.axhline(0.0, color="black", lw=0.8)
        ax.set_title(label, fontsize=9)
        ax.set_xlabel("x")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(r"$D(x)=F(x)-x$")
    fig.suptitle("M1A.1 — representative cumulative-excess profiles "
                 "(dotted: argmax boundary)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path)
    plt.close(fig)


# --------------------- M1B tolerance figures ---------------------


def plot_stopping_vs_epsilon(
    case_rows: list[dict],
    K_list: list[int],
    eps_entries: list[dict],
    *,
    path: Path,
) -> None:
    """Stopping round vs epsilon for each K (categorical x-axis, log-y stops)."""
    fig, ax = plt.subplots(figsize=(8.2, 5.0), dpi=150)
    colors = {50: "tab:green", 100: "tab:blue", 200: "tab:orange", 400: "tab:red"}
    xs = [e["epsilon"] for e in eps_entries]
    xpos = {e: i for i, e in enumerate(xs)}
    for K in sorted(K_list):
        pts = [(xpos[r["epsilon"]], r["stopping_round"]) for r in case_rows
               if r["K"] == K and r["stopping_round"] != "none"]
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=6,
                color=colors.get(K), label=f"K={K}")
    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels([f"{e:g}\n({e * 100:g}%)" for e in xs])
    ax.set_xlabel(r"tolerance $\epsilon$ (cumulative excess deadband)")
    ax.set_ylabel(r"stopping round $t_{\epsilon}$")
    ax.set_yscale("log")
    ax.set_title("M1B — finite tolerance removes the K-dependence of exact stopping\n"
                 r"($\epsilon=0$ row uses the canonical 1e-12 tolerance; T_max=5000)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_regret_vs_epsilon(
    case_rows: list[dict],
    K_list: list[int],
    eps_entries: list[dict],
    *,
    path: Path,
) -> None:
    """Uniformity regret at stopping vs epsilon, per K."""
    fig, ax = plt.subplots(figsize=(8.2, 5.0), dpi=150)
    colors = {50: "tab:green", 100: "tab:blue", 200: "tab:orange", 400: "tab:red"}
    xs = [e["epsilon"] for e in eps_entries]
    xpos = {e: i for i, e in enumerate(xs)}
    for K in sorted(K_list):
        pts = [(xpos[r["epsilon"]], float(r["uniformity_regret"]))
               for r in case_rows
               if r["K"] == K and r["stopping_round"] != "none"]
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=6,
                color=colors.get(K), label=f"K={K}")
    ax.axhline(0.0, color="black", lw=0.8, ls=":")
    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels([f"{e:g}\n({e * 100:g}%)" for e in xs])
    ax.set_xlabel(r"tolerance $\epsilon$")
    ax.set_ylabel(r"uniformity regret $=(U_{\rm stop}-U_{\rm best})/U_{\rm best}$")
    ax.set_yscale("symlog", linthresh=0.01)
    ax.set_title("M1B — how far the stopping state is from the same trajectory's\n"
                 "best achievable uniformity (0 = stops exactly at the U-optimum)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_m1b_reference(
    references: dict[int, dict],
    K_list: list[int],
    eps_entries: list[dict],
    *,
    path: Path,
) -> None:
    """Reference D_max trajectories with the epsilon deadband levels."""
    fig, ax = plt.subplots(figsize=(8.6, 5.2), dpi=150)
    colors = {50: "tab:green", 100: "tab:blue", 200: "tab:orange", 400: "tab:red"}
    for K in sorted(K_list):
        ref = references[K]
        stop = ref["stopped_at"] if ref["stopped_at"] is not None else len(ref["rows"]) - 1
        t = np.arange(stop + 1)
        ax.plot(t, ref["D_max_seq"][:stop + 1], lw=1.0, color=colors.get(K),
                label=f"K={K} (stops t={ref['stopped_at']})")
    for e in eps_entries:
        if e["epsilon"] > 0:
            ax.axhline(e["implemented_tolerance"], ls=":", lw=0.9,
                       label=fr"$\epsilon={e['epsilon']:g}$")
    ax.set_yscale("log")
    ax.set_xlabel("round t (pre-stop segment of the exact reference runs)")
    ax.set_ylabel(r"$D_{\max}(t)$ (log)")
    ax.set_title("M1B — first-passage view: stopping time = first crossing of the\n"
                 r"reference $D_{\max}(t)$ through $\epsilon$ (monotone in $\epsilon$ by construction)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------- M1B.1 metric-scale audit figures ---------------------


def plot_t_best_vs_K(
    opt_rows: list[dict],
    K_list: list[int],
    B_list: list[int],
    *,
    path: Path,
) -> None:
    """Reference-optimum timing vs simulation resolution, per observation scale."""
    fig, ax = plt.subplots(figsize=(8.2, 5.0), dpi=150)
    colors = {5: "tab:green", 10: "tab:orange", 25: "tab:red", 50: "tab:purple"}
    for B in B_list:
        pts = [(K, next(r["t_best"] for r in opt_rows if r["K"] == K and r["B"] == B))
               for K in K_list]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=6,
                color=colors.get(B), label=f"fixed coarse B={B}")
    fine = [(K, next(r["t_best"] for r in opt_rows if r["K"] == K and r["B"] == K))
            for K in K_list]
    ax.plot([p[0] for p in fine], [p[1] for p in fine], "s--", ms=7,
            color="black", label="fine grid (B=K)")
    ax.set_xlabel("simulation resolution K")
    ax.set_ylabel(r"reference optimum round $t_{\rm best}(K,B)$")
    ax.set_yscale("log")
    ax.set_title("M1B.1 — when is the optimum? fixed observation scales are\n"
                 "K-stable while the fine-grid optimum drifts with K")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_retention_vs_epsilon(
    stop_rows: list[dict],
    K_list: list[int],
    B_list: list[int],
    epsilons: list[float],
    *,
    path: Path,
    reference_K: int = 100,
) -> None:
    """Improvement retention at the M1B stopping state vs epsilon, per scale."""
    fig, ax = plt.subplots(figsize=(8.6, 5.2), dpi=150)
    colors = {5: "tab:green", 10: "tab:orange", 25: "tab:red", 50: "tab:purple",
              400: "tab:blue", 200: "tab:brown", 100: "black"}
    xs = range(len(epsilons))
    for B in B_list:
        pts = [next((float(r["improvement_retention"]) for r in stop_rows
                     if r["K"] == reference_K and r["epsilon"] == e and r["B"] == B),
                    np.nan) for e in epsilons]
        ax.plot(list(xs), pts, "o-", ms=6, color=colors.get(B),
                label=f"K={reference_K}, B={B}")
    # fine-grid line per K (thin)
    for K in K_list:
        pts = [next((float(r["improvement_retention"]) for r in stop_rows
                     if r["K"] == K and r["epsilon"] == e and r["B"] == K),
                    np.nan) for e in epsilons]
        ax.plot(list(xs), pts, "s--", ms=4, lw=1.0, alpha=0.65,
                color=colors.get(K), label=f"K={K}, B=K (fine)")
    for level, ls in [(0.99, ":"), (0.95, "--"), (0.90, "-.")]:
        ax.axhline(level, color="gray", ls=ls, lw=0.8)
        ax.annotate(f"{int(level * 100)}%", xy=(len(epsilons) - 0.9, level),
                    fontsize=7, color="gray", va="bottom")
    ax.set_xticks(list(xs))
    ax.set_xticklabels([f"{e:g}" for e in epsilons])
    ax.set_xlabel(r"tolerance $\epsilon$ (0 = exact M1A)")
    ax.set_ylabel("improvement retention at stopping state")
    ax.set_ylim(-0.02, 1.05)
    ax.set_title(f"M1B.1 — share of achievable uniformity improvement kept at stop\n"
                 f"(solid: coarse scales at K={reference_K}; dashed: fine grid per K; "
                 "reference lines are reading aids only)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_regret_vs_loss(
    stop_rows: list[dict],
    B_list: list[int],
    *,
    path: Path,
) -> None:
    """Residual regret vs improvement loss: the small-denominator contrast."""
    fig, ax = plt.subplots(figsize=(7.6, 5.4), dpi=150)
    scales = B_list + ["fine"]
    for B in scales:
        label = f"fine (B=K)" if B == "fine" else f"B={B}"
        pts = [(float(r["improvement_loss"]), float(r["residual_regret"]))
               for r in stop_rows
               if (r["scale_label"] == label
                   or (B == "fine" and r["scale_label"].startswith("fine")))]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=5,
                alpha=0.8, label=label)
    ax.set_xlabel("improvement loss (share of achievable gain thrown away)")
    ax.set_ylabel("residual regret (relative to U_best — small denominator!)")
    ax.set_yscale("log")
    ax.set_title("M1B.1 — same stopping states, two readings:\n"
                 "residual regret inflates as U_best approaches zero")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_coarse_profiles(
    profiles: dict[str, np.ndarray],
    *,
    B: int,
    path: Path,
) -> None:
    """Coarse macro-bin density profiles at a few representative states."""
    fig, axes = plt.subplots(1, len(profiles), figsize=(3.8 * len(profiles), 3.6),
                             dpi=150, sharey=True)
    if len(profiles) == 1:
        axes = [axes]
    xs = np.arange(1, B + 1) / B
    for ax, (label, P) in zip(axes, profiles.items()):
        ax.bar(xs, P, width=0.9 / B, color="tab:blue", alpha=0.7)
        ax.axhline(1.0 / B, color="black", lw=1.0, ls="--")
        ax.set_title(label, fontsize=9)
        ax.set_xlabel("macro-bin x (width 1/B)")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(f"macro-bin mass (B={B})")
    fig.suptitle("M1B.1 — coarse observation-scale density profiles "
                 "(dashed: perfectly uniform)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path)
    plt.close(fig)


# --------------------- M1B.2 micro-scaling diagnostic figures ---------------------


def plot_action_vs_time(runs: dict[int, dict], K_list: list[int], *, path: Path) -> None:
    """Gross action sizes vs normalized active time, per K (3 panels)."""
    colors = {50: "tab:green", 100: "tab:blue", 200: "tab:orange", 400: "tab:red"}
    fig, axes = plt.subplots(3, 1, figsize=(8.8, 8.6), dpi=150, sharex=True)
    quantities = [("moved_mass", "moved mass per round"),
                  ("net_export_mass", "net export per round"),
                  ("L1_state_step", r"$\|m_{t+1}-m_t\|_1$")]
    for ax, (key, label) in zip(axes, quantities):
        for K in sorted(K_list):
            diag = runs[K]["diag"]
            stop = runs[K]["stop"]
            frac = [d["t"] / stop for d in diag]
            ax.plot(frac, [d[key] for d in diag], lw=0.6, alpha=0.8,
                    color=colors.get(K), label=f"K={K}")
        medians = [np.median([d[key] for d in runs[K]["diag"]]) for K in sorted(K_list)]
        ax.axhline(np.mean(medians), color="black", ls=":", lw=1.0,
                   label=f"cross-K median ≈ {np.mean(medians):.3f}")
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, loc="upper right")
    axes[-1].set_xlabel("fraction of active trajectory (0 = start, 1 = exact stop)")
    axes[0].set_title("M1B.2 — gross action sizes do NOT shrink: stationary churn\n"
                      "at every resolution (micro-correction hypothesis falsified)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_churn_vs_D_max(
    runs: dict[int, dict], K_list: list[int], *, path: Path
) -> None:
    """Moved mass / net export vs current D_max: churn is independent of excess."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), dpi=150, sharex=True)
    colors = {50: "tab:green", 100: "tab:blue", 200: "tab:orange", 400: "tab:red"}
    for K in sorted(K_list):
        diag = runs[K]["diag"]
        Dm = [d["D_max"] for d in diag]
        axes[0].plot(Dm, [d["moved_mass"] for d in diag], "o", ms=2, alpha=0.5,
                     color=colors.get(K), label=f"K={K}")
        axes[1].plot(Dm, [d["net_export_mass"] for d in diag], "o", ms=2, alpha=0.5,
                     color=colors.get(K), label=f"K={K}")
    for ax, lab in zip(axes, ["moved mass", "net export mass"]):
        ax.set_xscale("log")
        ax.set_xlabel(r"$D_{\max}(t)$ (log)")
        ax.set_ylabel(f"{lab} per round")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle("M1B.2 — per-round action size vs current excess level:\n"
                 "no contraction — the churn continues at every excess level",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path)
    plt.close(fig)


def plot_reference_D_max_bands(
    runs: dict[int, dict],
    K_list: list[int],
    eps_list: list[float],
    *,
    path: Path,
) -> None:
    """Reference D_max trajectories (pre-stop) with epsilon levels and stops."""
    fig, ax = plt.subplots(figsize=(9.0, 5.2), dpi=150)
    colors = {50: "tab:green", 100: "tab:blue", 200: "tab:orange", 400: "tab:red"}
    for K in sorted(K_list):
        diag = runs[K]["diag"]
        ax.plot([d["t"] for d in diag], [d["D_max"] for d in diag], lw=0.7,
                color=colors.get(K), label=f"K={K} (stop t={runs[K]['stop']})")
        ax.axvline(runs[K]["stop"], color=colors.get(K), ls=":", lw=0.8)
    for eps in eps_list:
        ax.axhline(eps, color="gray", ls="--", lw=0.7)
        ax.annotate(fr"$\epsilon={eps:g}$", xy=(0.99, eps), xycoords=("axes fraction", "data"),
                    fontsize=7, va="bottom", ha="right", color="gray")
    ax.set_yscale("log")
    ax.set_xlabel("round t (active phase)")
    ax.set_ylabel(r"$D_{\max}(t)$ (log)")
    ax.set_title("M1B.2 — the exact controller fluctuates in a K-robust band;\n"
                 "exact stop = rare dip through zero, M1B stops at an ordinary dip")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_tail_vs_K(
    tail_rows: list[dict],
    over_rows: list[dict],
    *,
    path: Path,
) -> None:
    """Exact-stop tail length (after D_max thresholds) vs K, with M1B stops."""
    fig, ax = plt.subplots(figsize=(8.2, 5.0), dpi=150)
    Ks = [r["K"] for r in tail_rows]
    for thr in ["0.02", "0.01", "0.005"]:
        tails = [float(r[f"tail_over_K_{thr}"]) for r in tail_rows]
        ax.plot(Ks, tails, "o-", ms=6, label=f"tail after D_max≤{thr}, divided by K")
    for r in over_rows:
        if float(r["epsilon"]) > 0:
            ax.plot(r["K"], r["t_stop"] / r["K"], "^", ms=8, mfc="none", alpha=0.7,
                    color="gray")
    ax.set_xlabel("simulation resolution K")
    ax.set_ylabel("rounds (normalized by K)")
    ax.set_yscale("log")
    ax.set_title("M1B.2 — the exact-stopping tail is O(K) (micro-scale in wait time,\n"
                 "not in action size); triangles: M1B finite-epsilon stops (K-robust)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_a_vs_b_trajectories(runs: dict, K_primary: int, *, path) -> None:
    """M1C: Variant A (canonical q kernel) vs Variant B (deficit-fill kernel)
    at the primary resolution — uniformity, D_max, moved mass and boundary."""
    fig, axes = plt.subplots(4, 1, figsize=(8.8, 10.5), dpi=150, sharex=True)
    style = {"A": ("tab:blue", "Variant A — fixed biased q (canonical)"),
             "B": ("tab:red", "Variant B — target-aware deficit fill")}
    for variant, (color, label) in style.items():
        v = runs[(variant, K_primary)]
        stop = v["stop"]
        t = np.arange(stop)
        a_t = [r["a_t"] for r in v["run"]["rows"][:stop]]
        axes[0].plot(t, v["Ud"][:stop], color=color, lw=1.0, label=label)
        axes[1].plot(t, v["Dm"][:stop], color=color, lw=1.0)
        axes[2].plot(t, v["moved"][:stop], color=color, lw=1.0)
        axes[3].plot(t, a_t, color=color, lw=1.0, marker="o", ms=3)
    axes[0].set_yscale("log")
    axes[0].set_ylabel(r"$U_{\rm density} = K \cdot U$")
    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"$D_{\max}(t)$")
    axes[2].set_ylabel("moved mass per round")
    axes[3].set_ylabel(r"sweep boundary $a_t$")
    axes[3].set_xlabel("round t (active phase)")
    for ax in axes:
        ax.grid(True, alpha=0.3)
    axes[0].legend(fontsize=8, loc="upper right")
    axes[0].set_title(f"M1C correction-kernel ablation, K={K_primary}: selection rule "
                      "fixed, only the redistribution kernel changes")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_a_vs_b_summary(runs: dict, K_list: list, *, path) -> None:
    """M1C cross-K summary: exact stopping time and post-stop uniformity."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), dpi=150)
    style = {"A": ("tab:blue", "Variant A — fixed biased q (canonical)"),
             "B": ("tab:red", "Variant B — target-aware deficit fill")}
    for variant, (color, label) in style.items():
        Ks = sorted(K_list)
        stops = [runs[(variant, K)]["stop"] for K in Ks]
        u_stop = [runs[(variant, K)]["Ud"][runs[(variant, K)]["stop"] - 1]
                  for K in Ks]
        axes[0].plot(Ks, stops, "o-", color=color, label=label)
        axes[1].plot(Ks, u_stop, "o-", color=color, label=label)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("K")
    axes[0].set_ylabel("exact stopping round $t_{stop}$")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("K")
    axes[1].set_ylabel(r"$U_{\rm density}$ at the exact stop")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("M1C — churn lifetime collapses under the deficit-fill kernel\n"
                 "(same selection rule, same alpha, same stopping rules)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(path)
    plt.close(fig)


def plot_aub_trajectories(runs: dict, K_primary: int, *, path) -> None:
    """M1C.1: uniformity trajectories for A (near-biased q), U (fixed uniform
    q) and B (deficit-fill oracle) at the primary resolution."""
    fig, ax = plt.subplots(figsize=(8.8, 5.2), dpi=150)
    style = {"A": ("tab:blue", "Variant A — fixed near-biased $q$ (canonical)"),
             "U": ("tab:green", "Variant U — fixed uniform $q$ (target-matched, deficit-unaware)"),
             "B": ("tab:red", "Variant B — deficit-fill oracle")}
    for variant, (color, label) in style.items():
        v = runs[(variant, K_primary)]
        stop = v["stop"]
        ax.plot(np.arange(stop + 1), v["Ud"][:stop + 1], color=color, lw=1.2,
                marker="o", ms=3, label=label)
    ax.set_yscale("log")
    ax.set_xlabel("round t (states 0 .. exact stop)")
    ax.set_ylabel(r"$U_{\rm density} = K \cdot U$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title(f"M1C.1 fixed-uniform kernel ablation, K={K_primary}: "
                 "uniformity trajectories\n(selection rule fixed; only the "
                 "redistribution law differs)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_aub_kernel_windows(runs: dict, K_primary: int, *, path) -> None:
    """M1C.1: D_max and per-round moved mass for A / U / B at the primary
    resolution (A's churn window dwarfs U's 6-round and B's 3-round lives)."""
    fig, axes = plt.subplots(2, 1, figsize=(8.8, 7.6), dpi=150, sharex=True)
    style = {"A": ("tab:blue", "A — fixed near-biased $q$"),
             "U": ("tab:green", "U — fixed uniform $q$"),
             "B": ("tab:red", "B — deficit-fill oracle")}
    for variant, (color, label) in style.items():
        v = runs[(variant, K_primary)]
        stop = v["stop"]
        axes[0].plot(np.arange(stop), v["Dm"][:stop], color=color, lw=1.0,
                     label=label)
        axes[1].plot(np.arange(stop), v["moved"][:stop], color=color, lw=1.0)
    axes[0].set_yscale("log")
    axes[0].set_ylabel(r"$D_{\max}(t)$")
    axes[1].set_ylabel("moved mass per round")
    axes[1].set_xlabel("round t (active phase)")
    for ax in axes:
        ax.grid(True, alpha=0.3)
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle(f"M1C.1, K={K_primary}: removing the near-bias removes the "
                 "churn window\n(A churns for hundreds of rounds; U and B "
                 "terminate monotonically)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path)
    plt.close(fig)


def plot_R_threshold(res: dict, K_primary: int, *, path) -> None:
    """M1C.4: rebound-pressure ratio R_t vs normalized active time, with the
    R = 1 rebound threshold (K=100)."""
    fig, ax = plt.subplots(figsize=(9.2, 4.4), dpi=150)
    per = res[(K_primary, "per")]
    stop = res[(K_primary, "stop")]
    t = np.arange(stop) / stop
    R = np.array([float(r["R_t"]) for r in per])
    reb = np.array([r["rebound"] for r in per], dtype=bool)
    ax.plot(t, R, lw=0.5, color="tab:blue", label=r"$R_t$ (max reinjection / margin)")
    ax.plot(t[reb], R[reb], "o", ms=1.5, color="tab:red", label="rebound ($R_t > 1$)")
    ax.axhline(1.0, color="black", ls="--", lw=1.0, label=r"$R = 1$ threshold")
    onset = res[(K_primary, "onset")] / stop
    ax.axvline(onset, color="gray", ls=":", lw=1.0,
               label="churn-window onset (M1B.2 divider)")
    ax.set_xlabel(r"normalized active time $t / t_{stop}$")
    ax.set_ylabel(r"$R_t$")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title(f"M1C.4 — rebound-pressure ratio, canonical A, K={K_primary}: "
                 "rebound = the reinjection/margin race crossing 1")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_R_crossK(res: dict, K_list: list, *, path) -> None:
    """M1C.4: churn-window R_t quantile bands per K (cross-K stability)."""
    fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=150)
    colors = {50: "tab:green", 100: "tab:blue", 200: "tab:orange", 400: "tab:red"}
    for K in K_list:
        R = np.sort(res[(K, "R_w")])
        n = len(R)
        qs = [R[int(q * (n - 1))] for q in [0.05, 0.25, 0.5, 0.75, 0.95]]
        x = K * np.ones(5)
        ax.plot(x, qs, "o-", color=colors[K], label=f"K={K}")
    ax.set_xscale("log")
    ax.set_xlabel("K (p5 / p25 / p50 / p75 / p95 of churn-window $R_t$)")
    ax.set_ylabel(r"$R_t$ quantiles")
    ax.set_yscale("log")
    ax.axhline(1.0, color="black", ls="--", lw=1.0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title("M1C.4 — churn-window $R_t$ distribution across K:\n"
                 "the rebound fraction is P($R_t > 1$) for a stable band")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_phase_a_vs_R(res: dict, K_primary: int, *, path) -> None:
    """M1C.4: phase-style view — sweep boundary a_t vs R_t (K=100), colored
    by rebound flag."""
    fig, ax = plt.subplots(figsize=(7.6, 5.2), dpi=150)
    per = res[(K_primary, "per")]
    onset = res[(K_primary, "onset")]
    win = per[onset:]
    R = np.array([float(r["R_t"]) for r in win])
    a = np.array([float(r["a_t"]) for r in win])
    reb = np.array([r["rebound"] for r in win], dtype=bool)
    ax.plot(a[~reb], R[~reb], "o", ms=1.6, alpha=0.5, color="tab:blue",
            label=r"contraction rounds ($R_t \leq 1$)")
    ax.plot(a[reb], R[reb], "o", ms=1.6, alpha=0.5, color="tab:red",
            label="rebound rounds ($R_t > 1$)")
    ax.axhline(1.0, color="black", ls="--", lw=1.0)
    ax.set_xlabel(r"sweep boundary $a_t = j^*/K$")
    ax.set_ylabel(r"$R_t$")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title(f"M1C.4 — churn-window operating plane (a, R), K={K_primary}:\n"
                 "rebounds need both a strong-mismatch boundary and a thin margin")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_terminal_alignment(align_rows: list, *, path) -> None:
    """M1C.4: cross-K terminal alignment over tau = t_stop - t (last 40)."""
    fig, axes = plt.subplots(3, 1, figsize=(8.2, 8.4), dpi=150, sharex=True)
    colors = {50: "tab:green", 100: "tab:blue", 200: "tab:orange", 400: "tab:red"}
    for K in [50, 100, 200, 400]:
        rows = [r for r in align_rows if r["K"] == K]
        tau = np.array([r["tau"] for r in rows])
        Dm = np.array([float(r["D_max"]) for r in rows])
        R = np.array([float(r["R_t"]) for r in rows])
        M = np.array([float(r["moved_mass"]) for r in rows])
        axes[0].plot(tau, Dm, "o-", ms=3, lw=0.9, color=colors[K], label=f"K={K}")
        axes[1].plot(tau, R, "o-", ms=3, lw=0.9, color=colors[K])
        axes[2].plot(tau, M, "o-", ms=3, lw=0.9, color=colors[K])
    axes[0].set_yscale("log")
    axes[0].set_ylabel(r"$D_{\max}$")
    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"$R_t$")
    axes[1].axhline(1.0, color="black", ls="--", lw=0.8)
    axes[2].set_ylabel("moved mass $M$")
    axes[2].set_xlabel(r"$\tau = t_{stop} - t$ (rounds before the exact stop)")
    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.invert_xaxis()
    axes[0].legend(fontsize=8)
    fig.suptitle(r"M1C.4 — terminal alignment: the last rounds share the "
                 r"low-$D_{\max}$, small-prefix, sub-threshold-R structure",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path)
    plt.close(fig)


def plot_gate_scatter(per_rows: list, K_panels: list, *, path) -> None:
    """M1C.5: the boundary gate — a_t vs D_max with the analytic a_crit(D)
    curve; rebounds must lie strictly above the curve (necessary gate)."""
    fig, axes = plt.subplots(1, len(K_panels), figsize=(6.2 * len(K_panels), 5.2),
                             dpi=150, sharey=True)
    if len(K_panels) == 1:
        axes = [axes]
    D_line = np.linspace(1e-4, 0.25, 400)
    a_line = (1.0 - D_line + np.sqrt(D_line * (D_line + 2.0))) / 2.0
    for ax, K in zip(axes, K_panels):
        rows = [r for r in per_rows if int(r["K"]) == K and r["in_churn_window"]]
        a = np.array([float(r["a_t"]) for r in rows])
        D = np.array([float(r["D_max"]) for r in rows])
        reb = np.array([r["rebound"] for r in rows], dtype=bool)
        ax.plot(D[~reb], a[~reb], "o", ms=1.8, alpha=0.4, color="tab:blue",
                label="contraction")
        ax.plot(D[reb], a[reb], "o", ms=1.8, alpha=0.5, color="tab:red",
                label="rebound")
        ax.plot(D_line, a_line, "-", color="black", lw=1.4,
                label=r"$a_{\rm crit}(D)$ (necessary gate)")
        ax.set_xlabel(r"$D_{\max}(t)$")
        ax.set_title(f"K={K} (churn window)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="lower right")
    axes[0].set_ylabel(r"sweep boundary $a_t$")
    fig.suptitle("M1C.5 — boundary-gated rebound: every rebound lies above the "
                 "analytic necessary gate\n(necessary, not sufficient: blue "
                 "points above the curve are gate-open contractions)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(path)
    plt.close(fig)


def plot_gate_margin(per_rows: list, *, path) -> None:
    """M1C.5: distribution of the gate margin B_t = a_t - a_crit(D_max,t)
    for rebound vs contraction rounds (all K pooled, churn window)."""
    fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=150)
    rows = [r for r in per_rows if r["in_churn_window"]]
    B = np.array([float(r["B_t"]) for r in rows])
    reb = np.array([r["rebound"] for r in rows], dtype=bool)
    bins = np.linspace(min(B.min(), -0.1), max(B.max(), 0.1), 61)
    ax.hist(B[~reb], bins=bins, alpha=0.6, color="tab:blue",
            label=f"contraction (n={int((~reb).sum())})")
    ax.hist(B[reb], bins=bins, alpha=0.6, color="tab:red",
            label=f"rebound (n={int(reb.sum())})")
    ax.axvline(0.0, color="black", ls="--", lw=1.2, label="gate B = 0")
    ax.set_xlabel(r"gate margin $B_t = a_t - a_{\rm crit}(D_{\max,t})$")
    ax.set_ylabel("rounds")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title("M1C.5 — no rebounds below the gate (B <= 0); the B > 0 side\n"
                 "is mixed (necessary ≠ sufficient)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_profile_factor(per_rows: list, *, path) -> None:
    """M1C.5: profile factor P(x_r) = D(x_r)/D_max for rebound vs gate-open
    contraction rounds — the profile shape decides the mixed region."""
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=150)
    # restrict to the normal churn band: deep valleys (tiny D_max) make the
    # ratio P = D(x_r)/D_max explode without information
    rows = [r for r in per_rows
            if r["in_churn_window"] and r["P_profile"] != ""
            and float(r["D_max"]) > 0.005]
    P_r = np.array([float(r["P_profile"]) for r in rows if r["rebound"]])
    P_c = np.array([float(r["P_profile"]) for r in rows
                    if not r["rebound"] and r["gate_open"]])
    ax.hist(P_c, bins=40, alpha=0.6, color="tab:blue",
            label=f"gate-open contraction (n={len(P_c)})")
    ax.hist(P_r, bins=40, alpha=0.6, color="tab:red",
            label=f"rebound (n={len(P_r)})")
    ax.set_xlabel(r"profile factor $P(x_r) = D(x_r)/D_{\max}$ at the driving prefix")
    ax.set_ylabel("rounds")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title("M1C.5 — the gate opens the possibility; the D-profile shape\n"
                 "decides whether reinjection actually crosses the margin",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_gate_episode_alignment(res: dict, K_list: list, *, path) -> None:
    """M1C.6: B_t aligned at rebound rounds (t_r-2 .. t_r+4), median and
    quartiles per K — the gate exit / return cycle."""
    fig, axes = plt.subplots(1, len(K_list), figsize=(3.4 * len(K_list), 4.4),
                             dpi=150, sharey=True)
    colors = {50: "tab:green", 100: "tab:blue", 200: "tab:orange", 400: "tab:red"}
    for ax, K in zip(axes, K_list):
        ss = sorted(res[K]["align"])
        med = [float(np.median(res[K]["align"][s])) for s in ss]
        p25 = [float(np.percentile(res[K]["align"][s], 25)) for s in ss]
        p75 = [float(np.percentile(res[K]["align"][s], 75)) for s in ss]
        ax.fill_between(ss, p25, p75, alpha=0.3, color=colors[K])
        ax.plot(ss, med, "o-", color=colors[K], label=f"K={K}")
        ax.axhline(0.0, color="black", ls="--", lw=1.0)
        ax.axvline(0.0, color="gray", ls=":", lw=1.0)
        ax.set_xlabel("s (rounds from rebound)")
        ax.set_title(f"K={K}", fontsize=9)
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(r"$B = a - a_{\rm crit}(D_{\max})$")
    axes[0].legend(fontsize=8)
    fig.suptitle("M1C.6 — gate exit and return around rebounds:\n"
                 "B dips below 0 right after the rebound and re-crosses "
                 "within ~2 contraction rounds", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(path)
    plt.close(fig)


def plot_boundary_map(res: dict, K_primary: int, *, path) -> None:
    """M1C.6: empirical boundary map a_t -> a_{t+1}, rebound vs contraction."""
    fig, ax = plt.subplots(figsize=(6.4, 6.0), dpi=150)
    rows = res[K_primary]["map_rows"]
    a_t = np.array([float(r["a_t"]) for r in rows])
    a_n = np.array([float(r["a_next"]) for r in rows])
    reb = np.array([r["rebound"] for r in rows], dtype=bool)
    ax.plot(a_t[~reb], a_n[~reb], "o", ms=1.6, alpha=0.35, color="tab:blue",
            label="contraction rounds")
    ax.plot(a_t[reb], a_n[reb], "o", ms=1.6, alpha=0.5, color="tab:red",
            label="rebound rounds")
    ax.plot([0, 1], [0, 1], "-", color="black", lw=0.8, alpha=0.6)
    ax.set_xlabel(r"$a_t$ (sweep boundary at $t$)")
    ax.set_ylabel(r"$a_{t+1}$ (next argmax position)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title(f"M1C.6 — boundary map, K={K_primary}: large two-way jumps\n"
                 "with rebound rounds launching from above the gate",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_peak_gap_alignment(res: dict, K_primary: int, *, path) -> None:
    """M1C.6: competing-peak height gap H = D(j2) - D(j1) aligned at rebounds,
    plus contraction-round peak motion medians."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), dpi=150)
    K = K_primary
    stop = res[K]["stop"]
    peaks = res[K]["peaks"]
    H = np.array([ (p["D2"] - p["D1"]) if p["j2"] is not None else np.nan
                   for p in peaks ])
    reb = res[K]["reb"]
    t_r_list = res[K]["t_r_list"]
    ss = range(-3, 5)
    med, p25, p75 = [], [], []
    for s in ss:
        vals = [H[t_r + s] for t_r in t_r_list if 0 <= t_r + s < stop
                and not np.isnan(H[t_r + s])]
        med.append(np.median(vals))
        p25.append(np.percentile(vals, 25))
        p75.append(np.percentile(vals, 75))
    axes[0].fill_between(ss, p25, p75, alpha=0.3, color="tab:purple")
    axes[0].plot(ss, med, "o-", color="tab:purple")
    axes[0].axhline(0.0, color="black", ls="--", lw=1.0)
    axes[0].axvline(0.0, color="gray", ls=":", lw=1.0)
    axes[0].set_xlabel("s (rounds from rebound)")
    axes[0].set_ylabel(r"gap $H = D(j_2) - D(j_1)$")
    axes[0].set_title("competing-peak gap around rebounds", fontsize=9)
    axes[0].grid(True, alpha=0.3)
    # contraction peak motion
    d_act, d_comp = res[K]["d_act"], res[K]["d_comp"]
    axes[1].hist(d_act, bins=50, alpha=0.55, color="tab:blue",
                 label=rf"active peak $\Delta D_1$ (med {np.median(d_act):+.4f})")
    axes[1].hist(d_comp, bins=50, alpha=0.55, color="tab:green",
                 label=rf"competing peak $\Delta D_2$ (med {np.median(d_comp):+.4f})")
    axes[1].axvline(0.0, color="black", lw=0.8)
    axes[1].set_xlabel("one-round height change (contraction rounds)")
    axes[1].set_ylabel("rounds")
    axes[1].set_yscale("log")
    axes[1].legend(fontsize=8)
    axes[1].set_title("contraction lowers the active peak faster\n"
                      "than the competing peak", fontsize=9)
    for ax in axes:
        ax.grid(True, alpha=0.3)
    fig.suptitle(f"M1C.6 — peak competition, K={K}", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path)
    plt.close(fig)


def plot_terminal_vs_ordinary(res: dict, K_list: list, *, path) -> None:
    """M1C.6: B trajectory in the terminal episode vs ordinary episodes."""
    fig, ax = plt.subplots(figsize=(8.0, 4.6), dpi=150)
    colors = {50: "tab:green", 100: "tab:blue", 200: "tab:orange", 400: "tab:red"}
    for K in K_list:
        v = res[K]
        term = next(e for e in v["episodes"] if e["terminal"])
        s_range = range(0, term["tau_rebound"] + 1)
        Bt = [float(v["B"][term["t_r"] + s]) for s in s_range]
        ax.plot(s_range, Bt, "o-", ms=4, color=colors[K], label=f"K={K} terminal")
        ord_B = []
        eps = [e for e in v["episodes"] if not e["terminal"]]
        L = max(min(e["tau_rebound"] for e in eps), 3)
        for s in range(L + 1):
            vals = [float(v["B"][e["t_r"] + s]) for e in eps
                    if e["t_r"] + s < len(v["B"])]
            ord_B.append(np.median(vals))
        ax.plot(range(L + 1), ord_B, "--", color=colors[K], alpha=0.6,
                label=f"K={K} ordinary (median)")
    ax.axhline(0.0, color="black", ls="--", lw=1.0)
    ax.set_xlabel("s (rounds from episode-opening rebound)")
    ax.set_ylabel(r"$B = a - a_{\rm crit}(D_{\max})$")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    ax.set_title("M1C.6 — terminal vs ordinary episodes: the gate still "
                 "re-opens, but the trajectory\nenters the annihilable "
                 "geometry instead of the next rebound", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
