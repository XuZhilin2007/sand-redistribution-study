"""Rebound-witness spatial localization audit for M2 Gate 1
(analysis-only round; no new dynamics, no changes to Gate 1 raw results).

One command regenerates every audit CSV and the figure:

    python experiments/m2_gate1_generality/run_interpretation_audit.py

Research question (frozen): WHERE are the prefixes that actually trigger
rebound (witness) in KC / KS, and where is the closest-to-rebound risk
location in KW, which never crosses the threshold?

Witness definition (exact M1C.3 criterion quantities, per active round):

    E_t(j)   = M_t * Delta(j) - gamma_t(j)
    gamma_t(j) = D_max(t) - D_TM_new,t(j)
    j_witness(t) = argmax_j E_t(j)   (ties -> smallest j, np.argmax semantics)
    x_witness(t) = j_witness(t) / K

For KC / KS the witness is extracted on ACTUAL rebound rounds only
(max_j E_t(j) > 0); KW has zero rebounds, so its argmax is a RISK
location (max_j E_t(j) <= 0 on every round) and is never called a
witness.

Data provenance: the committed per-round CSV stores scalars only, and its
j_drive column is the argmax of the R_t RATIO (M1C.4 semantics), which is
NOT the difference-argmax witness. E_t(j) depends on the full D-profile,
so the frozen trajectories are regenerated deterministically from the
frozen config (identical run_m1a_deterministic calls, H(K)=10K) and
VERIFIED string-identical against the committed per-round diagnostics
(D_max_before, j_t, a_t, moved_mass, rebound, worst_excess, R_t) on every
active round BEFORE any witness derivation. Any mismatch aborts the audit.

Region classification (frozen): Near x < 1/3, Middle 1/3 <= x <= 2/3,
Far x > 2/3.
"""

from __future__ import annotations

import argparse
import csv
import json
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
    flat_top_distribution,
    weak_positive_distribution,
)

KERNELS = {"KC": None, "KW": weak_positive_distribution(),
           "KS": flat_top_distribution()}
K_LIST = [50, 100, 200, 400]
ALPHA, TOL = 0.25, 1e-12
AUDITED_ROUNDS = {"KC", "KS", "KW"}


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_csv_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def region_of(x: float) -> str:
    """Frozen three-way region split (Near / Middle / Far)."""
    if x < 1.0 / 3.0:
        return "Near"
    if x <= 2.0 / 3.0:
        return "Middle"
    return "Far"


def witness_argmax(E: np.ndarray) -> tuple[int, float]:
    """Deterministic witness: smallest-index argmax of E_t(j) (np.argmax
    semantics, identical to the repository's boundary tie-breaking)."""
    j = int(np.argmax(E))  # first (smallest) index on ties
    return j + 1, float(E[j])


def witness_argmax_prefix(E: np.ndarray, j_star: int) -> tuple[int, float]:
    """Supplementary prefix-restricted argmax (j <= j*): by the M1C.3
    suffix-impossibility remark any REBOUND-satisfying j must lie in the
    swept prefix, so this is the criterion-relevant risk location when the
    raw argmax lands on a structurally-impossible suffix position (where
    E(j) <= 0 always, and E(K) = -D_max exactly)."""
    j = int(np.argmax(E[:j_star]))
    return j + 1, float(E[j])


def E_profile(mass: np.ndarray, j_star: int, g_bins: np.ndarray) -> np.ndarray:
    """E_t(j) = M_t Delta(j) - gamma_t(j): the existing M1C.3/M1C.5
    `rebound_residual` helper (direct form) evaluated on the current state
    — reused verbatim, no new semantics."""
    return rebound_residual(mass, j_star, ALPHA, g_bins)


def quantiles(v: np.ndarray) -> dict:
    return {"count": int(len(v)),
            "mean": f"{float(v.mean()):.5f}" if len(v) else "",
            "median": f"{float(np.median(v)):.5f}" if len(v) else "",
            "p05": f"{float(np.percentile(v, 5)):.5f}" if len(v) else "",
            "p25": f"{float(np.percentile(v, 25)):.5f}" if len(v) else "",
            "p75": f"{float(np.percentile(v, 75)):.5f}" if len(v) else "",
            "p95": f"{float(np.percentile(v, 95)):.5f}" if len(v) else "",
            "min": f"{float(v.min()):.5f}" if len(v) else "",
            "max": f"{float(v.max()):.5f}" if len(v) else ""}


def region_shares(x_vals: list[float]) -> dict:
    n = len(x_vals)
    counts = {"Near": 0, "Middle": 0, "Far": 0}
    for x in x_vals:
        counts[region_of(x)] += 1
    return {"Near": f"{counts['Near'] / n:.5f}" if n else "",
            "Middle": f"{counts['Middle'] / n:.5f}" if n else "",
            "Far": f"{counts['Far'] / n:.5f}" if n else "",
            **{f"n_{k}": v for k, v in counts.items()}}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=Path(__file__).resolve().parent / "config.json",
    )
    args = parser.parse_args(argv)
    out_dir = args.config.resolve().parent / "results"
    q_near = DISTRIBUTIONS["near"]
    started = time.time()

    print("M2 Gate 1 interpretation audit — rebound-witness spatial localization")
    print("  E_t(j) = M_t Delta(j) - gamma_t(j); witness = smallest-index argmax")
    committed_rows = read_csv_rows(out_dir / "per_round_diagnostics.csv")
    committed = {(r["kernel"], int(r["K"]), int(r["t"])): r for r in committed_rows
                 if r["kernel"] in AUDITED_ROUNDS}

    rows_out: list[dict] = []
    move_rows: list[dict] = []
    x_by = {("KC", "rebound"): [], ("KS", "rebound"): [], ("KW", "risk"): []}
    R_by = {("KC", "rebound"): [], ("KS", "rebound"): [], ("KW", "risk"): []}
    dx_by = {"KC": [], "KS": []}

    # ---------------- regenerate + verify against committed reality ----------------
    print("\n=== regenerate frozen trajectories and verify vs committed per-round CSV ===")
    for kernel, respray in KERNELS.items():
        g = (respray if respray is not None else q_near).bin_probs(100)
        for K in K_LIST:
            g = (respray if respray is not None else q_near).bin_probs(K)
            run = run_m1a_deterministic(
                q_near, K=K, T_max=10 * K, alpha=ALPHA, stop_tolerance=TOL,
                redistribution="throw", redistribution_dist=respray)
            states = np.asarray(run["states"], dtype=float)
            rows = run["rows"]
            n = run["stopped_at"] if run["stopped_at"] is not None else 10 * K
            T = np.arange(1, K + 1) / K
            delta = np.cumsum(g) - T
            xw_series: list[float] = []      # per-round raw argmax-E locations
            xw_pref_series: list[float] = []  # per-round prefix-restricted argmax
            reb_idx_this_run: list[int] = []
            n_verified = 0
            for t in range(n):
                row = rows[t]
                c = committed[(kernel, K, t)]
                d_max = float(row["D_max"])
                j_star = int(row["j_t"])
                moved = float(row["removed_mass"])
                E = E_profile(states[t], j_star, g)
                j_w, max_E = witness_argmax(E)
                rebound_actual = float(rows[t + 1]["D_max"]) > d_max
                # criterion semantics reuse: max E > 0 iff actual rebound
                if (max_E > 0.0) != rebound_actual:
                    print(f"[FAIL] criterion mismatch {kernel} K={K} t={t}")
                    return 1
                # R_t recomputed with M1C.4 ratio semantics (verification only)
                gamma = d_max - target_matched_excess(states[t], j_star, ALPHA)
                pos = moved * np.maximum(delta, 0.0)
                ratio = np.where(gamma > 0.0, pos / np.where(gamma > 0.0, gamma, 1.0), 0.0)
                R = float(ratio.max())
                checks = {
                    "D_max_before": (f"{d_max:.10e}", c["D_max_before"]),
                    "j_t": (str(j_star), c["j_t"]),
                    "a_t": (f"{j_star / K:.5f}", c["a_t"]),
                    "moved_mass": (f"{moved:.10e}", c["moved_mass"]),
                    "rebound": (str(bool(rebound_actual)), c["rebound"]),
                    "worst_excess": (f"{max_E:.6e}", c["worst_excess"]),
                    "R_t": (f"{R:.10e}", c["R_t"]),
                }
                bad = [k for k, (got, exp) in checks.items() if got != exp]
                if bad:
                    print(f"[FAIL] regeneration mismatch {kernel} K={K} t={t}: "
                          f"{bad} (got {[(k, checks[k][0], checks[k][1]) for k in bad]})")
                    return 1
                n_verified += 1
                # witness/risk record (raw argmax = frozen definition) +
                # supplementary prefix-restricted argmax (criterion-relevant
                # when the raw argmax lands on a suffix position)
                j_pref, E_pref = witness_argmax_prefix(E, j_star)
                x_w = j_w / K
                dist_edge = min(abs(x_w - 1.0 / 3.0), abs(x_w - 2.0 / 3.0))
                kind = ("rebound" if rebound_actual else "contraction") \
                    if kernel != "KW" else "risk"
                rec = {
                    "kernel": kernel, "K": K, "t": t, "kind": kind,
                    "a_t": f"{j_star / K:.5f}",
                    "D_max": f"{d_max:.6f}",
                    "M_t": f"{moved:.8f}",
                    "j_witness": j_w,
                    "x_witness": f"{x_w:.5f}",
                    "x_minus_a": f"{x_w - j_star / K:.5f}",
                    "abs_x_minus_a": f"{abs(x_w - j_star / K):.5f}",
                    "max_E": f"{max_E:.8e}",
                    "R_t": f"{R:.8f}",
                    "region": region_of(x_w),
                    "argmax_on_suffix": bool(j_w > j_star),
                    "j_risk_prefix": j_pref,
                    "x_risk_prefix": f"{j_pref / K:.5f}",
                    "region_prefix": region_of(j_pref / K),
                    "dist_nearest_flattop_edge": f"{dist_edge:.5f}",
                }
                rows_out.append(rec)
                xw_series.append(x_w)
                xw_pref_series.append(j_pref / K)
                if rebound_actual:
                    x_by[(kernel, "rebound")].append(x_w)
                    R_by[(kernel, "rebound")].append(R)
                    move_rows.append({"kernel": kernel, "K": K, "t": t,
                                      "x_witness": f"{x_w:.5f}",
                                      "x_next": "", "delta_x": "",
                                      "abs_delta_x": ""})
                    reb_idx_this_run.append(len(move_rows) - 1)
                elif kernel == "KW":
                    x_by[("KW", "risk")].append(x_w)
                    R_by[("KW", "risk")].append(R)
            # resolve rebound->next movement within THIS (kernel, K) run
            for idx in reb_idx_this_run:
                mr = move_rows[idx]
                t0 = int(mr["t"])
                if t0 + 1 < len(xw_series):
                    mr["x_next"] = f"{xw_series[t0 + 1]:.5f}"
                    d = xw_series[t0 + 1] - float(mr["x_witness"])
                    mr["delta_x"] = f"{d:.5f}"
                    mr["abs_delta_x"] = f"{abs(d):.5f}"
                    if t0 + 1 < n:
                        dx_by[kernel].append(d)
            print(f"  [{kernel}] K={K}: {n_verified} active rounds verified "
                  f"string-identical (stop={run['stopped_at']})")

    # ---------------- summary tables ----------------
    summary_rows: list[dict] = []
    for kernel, sample in [("KC", "rebound"), ("KS", "rebound")]:
        for K in K_LIST:
            sel = [r for r in rows_out
                   if r["kernel"] == kernel and int(r["K"]) == K
                   and r["kind"] == "rebound"]
            x = np.array([float(r["x_witness"]) for r in sel])
            axa = np.array([float(r["abs_x_minus_a"]) for r in sel])
            summary_rows.append({
                "kernel": kernel, "K": K, "sample": "actual rebound witnesses",
                **quantiles(x),
                **region_shares(list(x)),
                "median_abs_x_minus_a": f"{float(np.median(axa)):.5f}" if len(axa) else "",
            })
    # KW risk: three frozen sample definitions, raw argmax AND
    # prefix-restricted argmax (criterion-relevant supplement)
    for K in K_LIST:
        sel = [r for r in rows_out if r["kernel"] == "KW" and int(r["K"]) == K]
        R = np.array([float(r["R_t"]) for r in sel])
        top25 = [r for r, rv in zip(sel, R) if rv >= float(np.percentile(R, 75))]
        closest = [max(sel, key=lambda r: float(r["max_E"]))]
        for loc_field, tag_prefix in [("x_witness", "raw argmax"),
                                      ("x_risk_prefix", "prefix-restricted")]:
            for tag, subset in [("all active rounds (risk)", sel),
                                ("top 25% R_t rounds (risk)", top25),
                                ("closest-to-rebound round (risk)", closest)]:
                x = np.array([float(r[loc_field]) for r in subset])
                axa = np.abs(x - np.array([float(r["a_t"]) for r in subset]))
                summary_rows.append({
                    "kernel": "KW", "K": K,
                    "sample": f"{tag} [{tag_prefix}]", **quantiles(x),
                    **region_shares(list(x)),
                    "median_abs_x_minus_a": f"{float(np.median(axa)):.5f}" if len(axa) else "",
                })
    write_csv(out_dir / "witness_localization_rows.csv", list(rows_out[0]), rows_out)
    write_csv(out_dir / "witness_localization_summary.csv",
              list(summary_rows[0]), summary_rows)
    write_csv(out_dir / "witness_next_movement.csv",
              list(move_rows[0]), move_rows)

    # ---------------- figure 1: location distributions ----------------
    fig, axes = plt.subplots(3, 1, figsize=(9, 7.5), sharex=True)
    panels = [
        ("KC rebound witnesses (pooled 4 K, n=%d)", x_by[("KC", "rebound")], "#1f77b4"),
        ("KS rebound witnesses (pooled 4 K, n=%d)", x_by[("KS", "rebound")], "#2ca02c"),
        ("KW risk locations, all active rounds (pooled 4 K, n=%d)", x_by[("KW", "risk")], "#ff7f0e"),
    ]
    bins = np.linspace(0.0, 1.0, 21)
    for ax, (label, vals, color) in zip(axes, panels):
        ax.hist(vals, bins=bins, color=color, edgecolor="black", lw=0.4)
        ax.set_ylabel("count")
        ax.set_title(label % len(vals), fontsize=10)
        for b in (1 / 3, 2 / 3):
            ax.axvline(b, color="gray", ls="--", lw=0.8)
    axes[-1].set_xlabel("x = j/K  (Near < 1/3, Middle, Far > 2/3)")
    fig.suptitle("M2 Gate 1 interpretation audit — witness / risk prefix locations")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_dir / "witness_location_distribution.png", dpi=150)
    plt.close(fig)

    # ---------------- figure 2: rebound -> next-round movement ----------------
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2), sharex=True, sharey=True)
    for ax, kernel, color in zip(axes, ["KC", "KS"],
                                 ["#1f77b4", "#2ca02c"]):
        xs = [float(r["x_witness"]) for r in move_rows if r["kernel"] == kernel
              and r["x_next"] != ""]
        xn = [float(r["x_next"]) for r in move_rows if r["kernel"] == kernel
              and r["x_next"] != ""]
        ax.scatter(xs, xn, s=9, alpha=0.5, color=color)
        ax.plot([0, 1], [0, 1], color="gray", lw=0.8, ls="--")
        ax.set_title(f"{kernel}: x_witness(t) vs x_next(t+1), n={len(xs)}",
                     fontsize=10)
        ax.set_xlabel("x_witness(t)")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    axes[0].set_ylabel("x_next(t+1)")
    fig.suptitle("Rebound-to-next-round witness/risk relocation")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_dir / "witness_next_movement.png", dpi=150)
    plt.close(fig)

    # ---------------- console summary ----------------
    print("\n=== §9 comparison table (actual rebound witnesses for KC/KS; "
          "KW = risk locations, NOT witnesses) ===")
    hdr = (f"{'kern':>4} {'K':>4} {'n':>5} {'x med':>7} {'p05-p95':>17} "
           f"{'Near':>7} {'Mid':>7} {'Far':>7} {'|x-a| med':>9}")
    print(hdr)
    for r in summary_rows:
        if "all active" in r["sample"] or r["sample"].startswith("actual"):
            print(f"{r['kernel']:>4} {r['K']:>4} {r['count']:>5} {r['median']:>7} "
                  f"{r['p05']}-{r['p95']:>11} {r['Near']:>7} {r['Middle']:>7} "
                  f"{r['Far']:>7} {r['median_abs_x_minus_a']:>9}")
    print("\nKW per-K closest-to-rebound rounds (raw argmax; prefix-restricted "
          "in brackets):")
    for K in K_LIST:
        sel = [r for r in rows_out if r["kernel"] == "KW" and int(r["K"]) == K]
        c = max(sel, key=lambda r: float(r["max_E"]))
        print(f"  K={K}: t={c['t']} x_risk={c['x_witness']} "
              f"(prefix {c['x_risk_prefix']}) a_t={c['a_t']} R_t={c['R_t']} "
              f"-gap(max E)={float(c['max_E']):.3e} region={c['region']} "
              f"region_prefix={c['region_prefix']} "
              f"argmax_on_suffix={c['argmax_on_suffix']}")
    n_sfx = sum(1 for r in rows_out if r["kernel"] == "KW" and r["argmax_on_suffix"])
    n_kw = sum(1 for r in rows_out if r["kernel"] == "KW")
    n_sfx_reb = sum(1 for r in rows_out
                    if r["kernel"] in ("KC", "KS") and r["kind"] == "rebound"
                    and r["argmax_on_suffix"])
    n_reb = sum(1 for r in rows_out if r["kind"] == "rebound")
    print(f"suffix-trivial raw argmax: KW {n_sfx}/{n_kw} rounds; "
          f"KC/KS rebound witnesses {n_sfx_reb}/{n_reb}")
    for kernel in ["KC", "KS"]:
        dx = np.array(dx_by[kernel])
        if len(dx):
            print(f"{kernel} rebound->next displacement: n={len(dx)} "
                  f"median |dx|={float(np.median(np.abs(dx))):.3f} "
                  f"p05/p95=[{float(np.percentile(dx, 5)):.3f},{float(np.percentile(dx, 95)):.3f}]")
    print(f"\naudit outputs written to: {out_dir}")
    print(f"runtime: {round(time.time() - started, 2)} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
