"""Independent exact audit of the uniform W0 proof and all 5,535 sweep points.

No experiment generator or symbolic kernel is imported.  Infinite-K coverage
comes from the analytic L1/variance proof, not from finite sampling.
"""
from collections import Counter
from fractions import Fraction as Q
import argparse
import hashlib
import json
from pathlib import Path
import time

from audit_v04_interval_v3 import (
    COMMON, AuditFailure, add, sub, mul, scale, total, power, trim, value,
    bernstein_coefficients, require, xo_cdf,
)


ROOT = Path(__file__).resolve().parents[1]
ALPHAS = (COMMON[0], Q(1, 4), COMMON[1])
SWEEP_KS = tuple(range(204, 2049))
REVIEW = ROOT / "experiments/v04_all_k/review_2026-10-23"


def loss_integral(pieces):
    result = ()
    for lo, hi, slope, intercept in pieces:
        d = sub(intercept, (Q(1),))
        result = add(result, total((
            scale(mul(mul(slope, slope), sub(power(hi, 3), power(lo, 3))), Q(1, 3)),
            mul(mul(slope, d), sub(power(hi, 2), power(lo, 2))),
            mul(mul(d, d), sub(hi, lo)),
        )))
    return result


def continuum_polynomials():
    alpha, rho = (Q(0), Q(1)), (Q(1), Q(-1))
    lam0 = scale(alpha, Q(3, 4))
    c = lambda x: (Q(x),)
    common = [(c(0), c(Q(1, 3)), Q(1, 2), True),
              (c(Q(1, 3)), c(Q(1, 2)), Q(3, 2), True),
              (c(Q(1, 2)), c(Q(2, 3)), Q(3, 2), False),
              (c(Q(2, 3)), c(Q(5, 6)), Q(2), False),
              (c(Q(5, 6)), c(1), Q(0), False)]
    first = []
    for lo, hi, gd, swept in common:
        factor = rho if swept else c(1)
        first.append((lo, hi, scale(factor, -2), add(scale(factor, 2), scale(lam0, gd))))
    answer = {"U0": c(Q(1, 3)), "U1": loss_integral(first)}
    for upper in (False, True):
        y = add(c(Q(1, 2)), scale(alpha, Q(3, 4) if upper else Q(9, 16)))
        gy = sub(scale(y, 2), c(Q(2, 3))) if upper else sub(scale(y, Q(3, 2)), c(Q(1, 3)))
        fs = sub(sub(scale(y, 2), mul(y, y)), mul(lam0, sub(c(1), gy)))
        lam1 = mul(alpha, fs)
        regions = [(c(0), c(Q(1, 3)), Q(1, 2), 2),
                   (c(Q(1, 3)), c(Q(1, 2)), Q(3, 2), 2)]
        if upper:
            regions += [(c(Q(1, 2)), c(Q(2, 3)), Q(3, 2), 1),
                        (c(Q(2, 3)), y, Q(2), 1), (y, c(Q(5, 6)), Q(2), 0)]
        else:
            regions += [(c(Q(1, 2)), y, Q(3, 2), 1),
                        (y, c(Q(2, 3)), Q(3, 2), 0),
                        (c(Q(2, 3)), c(Q(5, 6)), Q(2), 0)]
        regions += [(c(Q(5, 6)), c(1), Q(0), 0)]
        second = []
        for lo, hi, gd, sweeps in regions:
            original_factor = power(rho, sweeps)
            spray_factor = rho if sweeps else c(1)
            spray = add(mul(lam0, spray_factor), lam1)
            second.append((lo, hi, scale(original_factor, -2),
                           add(scale(original_factor, 2), scale(spray, gd))))
        answer["U2_upper" if upper else "U2_middle"] = loss_integral(second)
    return answer


def rounded_selector(k, alpha):
    u0 = Q(k // 2, k)
    fs0 = 2 * u0 - u0 * u0
    v = k * (Q(1, 2) + 3 * alpha * fs0 / 4)
    lower = v.numerator // v.denominator
    return lower if v - lower <= Q(1, 2) else lower + 1


def moment_replay(k, alpha):
    """Exact bin sums, using the proved selector on J, with O(1) regions."""
    s0, s1 = k // 2, rounded_selector(k, alpha)
    u0, u1, rho, h = Q(s0, k), Q(s1, k), 1 - alpha, Q(1, k)
    lam0 = alpha * (2 * u0 - u0 * u0)
    fs1 = 2 * u1 - u1 * u1 - lam0 * (1 - xo_cdf(u1))
    lam1 = alpha * fs1
    edges = {1, k + 1, s0 + 1, s1 + 1}
    for knot in (Q(1, 3), Q(2, 3), Q(5, 6)):
        j = (k * knot.numerator + knot.denominator - 1) // knot.denominator
        edges.update((j, j + 1))
    edges = sorted(x for x in edges if 1 <= x <= k + 1)
    losses = []
    for t in (1, 2):
        squared_sum = Q(0)
        for start, end in zip(edges[:-1], edges[1:]):
            last, n = end - 1, end - start
            gmass = xo_cdf(Q(start, k)) - xo_cdf(Q(start - 1, k))
            if t == 1:
                factor, beta = (rho if start <= s0 else Q(1)), lam0
            else:
                factor = rho ** 2 if start <= s0 else rho if start <= s1 else Q(1)
                beta = rho * lam0 + lam1 if start <= s1 else lam0 + lam1
            aa = factor * Q(2 * k + 1, k * k) + beta * gmass - h
            bb = -2 * factor / (k * k)
            si = Q((start + last) * n, 2)
            si2 = Q(last * (last + 1) * (2 * last + 1) -
                    (start - 1) * start * (2 * start - 1), 6)
            squared_sum += aa * aa * n + 2 * aa * bb * si + bb * bb * si2
        losses.append(k * squared_sum)
    return (Q(k * k - 1, 3 * k * k), *losses), (s0, s1)


def full_replay(k, alpha):
    """Independent full-state argmax replay, used only as a crosscheck."""
    p = [Q(2 * k - 2 * i + 1, k * k) for i in range(1, k + 1)]
    weights = [xo_cdf(Q(i, k)) - xo_cdf(Q(i - 1, k)) for i in range(1, k + 1)]
    losses, words, second_ties = [], [], None
    for t in range(3):
        losses.append(k * sum(((x - Q(1, k)) ** 2 for x in p), Q(0)))
        if t == 2:
            break
        prefix, discrepancies = Q(0), []
        for i, x in enumerate(p, 1):
            prefix += x
            discrepancies.append(prefix - Q(i, k))
        maximum = max(discrepancies)
        winners = [i for i, x in enumerate(discrepancies, 1) if x == maximum]
        require(maximum > Q(1, 10 ** 12), "crosscheck encountered stopping")
        selected = winners[0]
        if t == 1:
            second_ties = winners
        words.append(selected)
        mass = alpha * sum(p[:selected])
        p = [((1 - alpha) * x if i < selected else x) + mass * weights[i]
             for i, x in enumerate(p)]
    return tuple(losses), tuple(words), second_ties


def check_uniform_constants(polys):
    amax, hmax = Q(251, 1000), Q(1, 33)
    require(COMMON[1] < amax < Q(2, 7), "density/alpha bound")
    require(Q(34, 33) + 3 * amax / 2 < 2, "f1 density bound")
    kappa = Q(1, 2) + 3 * amax / 8 + 3 * amax * hmax / 16
    variance = lambda beta: (Q(1, 3) + 7 * beta / 4) * hmax + 21 * beta * beta / 16
    c1 = 4 * amax + variance(3 * amax / 4)
    c2 = 4 * amax + 8 * amax * kappa + variance(7 * amax / 4)
    require(c1 < Q(11, 10) and c2 < Q(5, 2), "uniform error constants")
    bounds = {"alpha_upper": amax, "h_upper": hmax, "kappa": kappa,
              "C1_raw": c1, "C1": Q(11, 10), "C2_raw": c2, "C2": Q(5, 2),
              "C1_plus_C2": Q(18, 5),
              "margin12_at_N145": Q(1, 40) - Q(18, 5 * 145),
              "margin02_at_N145": Q(9, 50) - Q(5, 2 * 145) - Q(1, 3 * 145 * 145)}
    require(bounds["margin12_at_N145"] > 0 and bounds["margin02_at_N145"] > 0,
            "large K margin")
    obligations = []
    for name, p in (("U0-U2M-9/50", sub(sub(polys["U0"], polys["U2_middle"]), (Q(9, 50),))),
                    ("U1-U2M-1/40", sub(sub(polys["U1"], polys["U2_middle"]), (Q(1, 40),)))):
        coeffs = bernstein_coefficients(p, *COMMON)
        require(all(x > 0 for x in coeffs), name)
        obligations.append({"name": name, "interval": [str(x) for x in COMMON],
                            "bernstein_coefficients": [str(x) for x in coeffs]})
    return {k: str(v) for k, v in bounds.items()}, obligations


def check_sweep(polys):
    path = ROOT / "experiments/v04_all_k/w0_sweep.jsonl"
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    require(Counter(r["K"] for r in rows) == Counter(SWEEP_KS), "sweep exact K coverage")
    minima, max_abs_r, max_at = {}, Q(-1), None
    count = 0
    for row in rows:
        require(set(row) == {"K", *(str(x) for x in ALPHAS)}, "sweep fields/alphas")
        k = row["K"]
        for alpha in ALPHAS:
            uu, words = moment_replay(k, alpha)
            expected = {"d02": uu[0] - uu[2], "d12": uu[1] - uu[2],
                        "R": k * (uu[2] - value(polys["U2_middle"], alpha))}
            cell = row[str(alpha)]
            require(set(cell) == {"d02", "d12", "R", "ok"} and cell["ok"] is True,
                    f"K={k}: sweep record structure/flag")
            for field, v in expected.items():
                require(isinstance(cell[field], str) and Q(cell[field]) == v,
                        f"K={k} alpha={alpha}: {field} mismatch")
            require(expected["d02"] > 0 and expected["d12"] > 0, "nonpositive W0 scan")
            for field in ("d02", "d12"):
                if field not in minima or expected[field] < minima[field][0]:
                    minima[field] = (expected[field], k, str(alpha))
            if abs(expected["R"]) > max_abs_r:
                max_abs_r, max_at = abs(expected["R"]), (k, str(alpha))
            count += 1
        if (k - 203) % 400 == 0:
            print(f"Sweep independently replayed {k - 203}/{len(SWEEP_KS)} K", flush=True)
    require(max_abs_r < Q(8681, 50000), "decimal upper bound 0.17362")
    return {"K_count": len(rows), "point_count": count,
            "minima": {key: {"value": str(v), "K": k, "alpha": a}
                       for key, (v, k, a) in minima.items()},
            "max_abs_R": str(max_abs_r), "max_abs_R_at": max_at,
            "strict_decimal_upper_bound": "0.17362"}


def check_small_complement():
    path = REVIEW / "finite_common_audit.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    require(report["status"] == "VERIFIED" and not report["failures"], "small-K audit status")
    require(set(range(6, 145)) <= set(report["verified_K"]), "small-K complement coverage")
    require(report["common_interval"] == [str(x) for x in COMMON], "small-K interval")
    # The generator/kernel hashes in the old audit are generation provenance;
    # the independent auditor never imports them.  Bind the actual mathematical
    # dependencies, allowing unrelated generator repairs after that audit.
    dependencies = ("tools/audit_v04_interval_v3.py",
                    "experiments/v04_w_interval/interval_certs_v3.jsonl",
                    "docs/XO_ALPHA_EXTENSION_PROOF.md",
                    "docs/S3_XO_BOUNDARY_THEOREM_PROOF.md")
    for p in dependencies:
        recorded_hashes = {name.replace("\\", "/"): value
                           for name, value in report["source_sha256"].items()}
        digest = recorded_hashes[p]
        require(hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == digest,
                f"small-K audit mathematical dependency changed: {p}")
    return {"from": 6, "through": 144, "dependency": str(path.relative_to(ROOT)),
            "checked_mathematical_dependencies": list(dependencies)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path,
                        help="write a new audit report; default only prints results")
    args = parser.parse_args()
    started = time.monotonic()
    polys = continuum_polynomials()
    source = ROOT / "experiments/v04_all_k/continuum_U_polys.json"
    recorded = {k: trim(Q(x) for x in v) for k, v in json.loads(source.read_text()).items()}
    require(polys == recorded, "continuum polynomial coefficient mismatch")
    constants, obligations = check_uniform_constants(polys)
    broad_obligations = []
    for branch, lo, hi in (("middle", Q(1, 5), Q(16, 63)),
                           ("upper", Q(16, 63), Q(7, 24))):
        for ref in ("U0", "U1"):
            coeffs = bernstein_coefficients(sub(polys[ref], polys["U2_" + branch]), lo, hi)
            require(all(x > 0 for x in coeffs), f"broad branchwise {branch} {ref}-U2")
            broad_obligations.append({"branch": branch, "name": ref + "-U2",
                                      "interval": [str(lo), str(hi)],
                                      "bernstein_coefficients": [str(x) for x in coeffs]})
    checks = []
    for k in (33, 53, 96, 145, 210, 238, 1000, 2048):
        for alpha in ALPHAS:
            uu, word = moment_replay(k, alpha)
            direct, actual_word, _ = full_replay(k, alpha)
            require(uu == direct and word == actual_word, f"full/moment replay K={k}, alpha={alpha}")
            checks.append({"K": k, "alpha": str(alpha), "word": word})
    small = check_small_complement()
    sweep = check_sweep(polys)
    at = Q(16, 63)
    switch = {"alpha": str(at), "middle_U2": str(value(polys["U2_middle"], at)),
              "upper_U2": str(value(polys["U2_upper"], at)),
              "upper_minus_middle": str(value(sub(polys["U2_upper"], polys["U2_middle"]), at)),
              "replays": []}
    for k in (42, 46, 84, 88, 420, 424, 4200, 4204):
        uu, word, ties = full_replay(k, at)
        expected = [27 * (k // 42), 29 * (k // 42)] if k % 42 == 0 else [29 * (k // 42) + 3]
        require(ties == expected, "switch-point subsequence branch")
        switch["replays"].append({"K": k, "word": word, "s1_maximizers": ties, "U2": str(uu[2])})
    sources = [Path(__file__).resolve(), ROOT / "tools/audit_v04_interval_v3.py",
               ROOT / "docs/v04_research/W0_UNIFORM_BOUND_REVIEW.md",
               ROOT / "docs/v04_research/LARGE_K_PREFIX_LEMMA_REVIEW.md",
               ROOT / "experiments/v04_all_k/w0_sweep.py",
               ROOT / "experiments/v04_all_k/w0_sweep.jsonl",
               ROOT / "experiments/v04_all_k/continuum_limit.py", source,
               REVIEW / "finite_common_audit.json"]
    report = {"status": "VERIFIED", "statement": "W0 for ALL K>=6 and ALL alpha in J",
              "J": [str(x) for x in COMMON], "uniform_error_domain_K_min": 33,
              "analytic_large_K_threshold": 145, "small_complement": small,
              "constants": constants, "bernstein_obligations": obligations,
              "broad_branchwise_W0_obligations": broad_obligations,
              "continuum_polynomials": {k: [str(x) for x in v] for k, v in polys.items()},
              "full_vector_crosschecks": checks, "sweep": sweep,
              "broad_window_switch_caveat": switch,
              "source_sha256": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sources},
              "analytic_dependencies": ["proved global prefix selector on J for K>=33",
                                        "exact bin-integral lift of two updates",
                                        "L1 contraction of fixed-cutoff positive Markov operator",
                                        "cutoff displacement bounds and density <=2",
                                        "exact cell variance with at most three kernel jumps",
                                        "Bernstein positivity on J and finite-complement interval certificates"],
              "elapsed_seconds": round(time.monotonic() - started, 3)}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n",
                               encoding="utf-8")
    print(f"VERIFIED: uniform W0 for all K>=6, all alpha in J; "
          f"{sweep['point_count']} exact sweep records independently replayed; "
          f"{len(checks)} full-vector crosschecks; {report['elapsed_seconds']:.2f}s", flush=True)


if __name__ == "__main__":
    main()
