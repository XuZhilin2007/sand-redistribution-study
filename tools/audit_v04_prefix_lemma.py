"""Exact constants and finite crosschecks for LARGE_K_PREFIX_LEMMA_REVIEW.md.

The infinite-K proof is the analytic argument in the document, not the finite
sample.  The only imported local module is the independent review auditor.
"""
from fractions import Fraction as Q
import argparse
import hashlib
import json
from pathlib import Path

from audit_v04_interval_v3 import COMMON, xo_cdf


ROOT = Path(__file__).resolve().parents[1]


def winners_after_one(k, alpha):
    s0 = k // 2
    p0 = [Q(2 * k - 2 * i + 1, k * k) for i in range(1, k + 1)]
    fs0 = sum(p0[:s0])
    lam = alpha * fs0
    prefix, maximum, winners = Q(0), None, []
    for i, p in enumerate(p0, 1):
        g = xo_cdf(Q(i, k)) - xo_cdf(Q(i - 1, k))
        prefix += ((1 - alpha) * p if i <= s0 else p) + lam * g
        discrepancy = prefix - Q(i, k)
        if maximum is None or discrepancy > maximum:
            maximum, winners = discrepancy, [i]
        elif discrepancy == maximum:
            winners.append(i)
    return winners, fs0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path,
                        help="write a new audit report; default only prints results")
    args = parser.parse_args()
    amin, amax = COMMON
    lam_min, lam_max = amin * Q(800, 1089), amax * Q(3, 4)
    peak_gap = lam_max * (4 - 21 * lam_max) / 48
    bounds = {
        "F0_s0_min": Q(800, 1089), "lambda_min": lam_min,
        "lambda_lower_coarse": Q(9, 50), "lambda_max": lam_max,
        "peak_gap_min": peak_gap, "rounding_loss_max": Q(1, 4 * 33 ** 2),
        "strict_global_gap": peak_gap - Q(1, 4 * 33 ** 2),
        "middle_grid_below_2over3": Q(2, 3) - (Q(1, 2) + 3 * lam_max / 4 + Q(1, 66)),
        "middle_grid_above_s0": 3 * Q(9, 50) / 4 - Q(1, 66),
        "prefix_low_derivative_min": Q(1, 3) - 4 * amax / 3 + Q(9, 100),
        "prefix_middle_derivative_factor": Q(3, 2) * Q(800, 1089) - 1,
        "middle_over_prefix_gap_min": (3 * Q(9, 50) / 4) ** 2 - Q(1, 4 * 33 ** 2),
    }
    if not (lam_min >= Q(9, 50) and Q(9, 50) > Q(2, 21) and lam_max < Q(4, 21)):
        raise AssertionError("lambda domain/monotonicity bounds")
    if not all(x > 0 for x in bounds.values()):
        raise AssertionError("nonpositive analytic bound")
    ks = list(range(33, 204)) + [224, 288, 300, 600, 999, 1000, 1001, 4095, 4096, 4097]
    ties = []
    for k in ks:
        for alpha in (amin, Q(1, 4), amax):
            actual, fs0 = winners_after_one(k, alpha)
            v = k * (Q(1, 2) + 3 * alpha * fs0 / 4)
            lower = v.numerator // v.denominator
            expected = ([lower, lower + 1] if v - lower == Q(1, 2)
                        else [lower if v - lower < Q(1, 2) else lower + 1])
            if actual != expected:
                raise AssertionError((k, alpha, expected, actual))
            if len(actual) > 1:
                ties.append({"K": k, "alpha": str(alpha), "winners": actual})
    small_ties = {k: actual for k in range(6, 33)
                  if len(actual := winners_after_one(k, Q(1, 4))[0]) > 1}
    if small_ties != {16: [10, 11], 32: [20, 21, 22]}:
        raise AssertionError(small_ties)
    sources = [Path(__file__).resolve(), ROOT / "tools/audit_v04_interval_v3.py",
               ROOT / "docs/v04_research/LARGE_K_PREFIX_LEMMA_REVIEW.md"]
    report = {
        "status": "EXACT_BOUND_CHECKS_AND_FINITE_CROSSCHECK_PASS",
        "infinite_K_dependency": "analytic proof in LARGE_K_PREFIX_LEMMA_REVIEW.md, sections 2-5",
        "source_sha256": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sources},
        "analytic_bounds": {k: str(v) for k, v in bounds.items()},
        "large_K_crosscheck_K": ks, "large_K_crosscheck_alphas": [str(amin), "1/4", str(amax)],
        "large_K_crosscheck_pairs": len(ks) * 3, "ties_seen": ties,
        "small_K_crosscheck_K": list(range(6, 33)), "small_K_ties": small_ties,
        "total_crosscheck_pairs": len(ks) * 3 + 27,
    }
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: exact analytic bounds, {report['total_crosscheck_pairs']} finite crosschecks; "
          f"small ties={small_ties}")


if __name__ == "__main__":
    main()
