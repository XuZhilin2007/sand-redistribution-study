"""Exact uniform tail proof for K>=204, alpha in J.

The proof identities and analytic dependencies are documented in
UNIFORM_TAIL_CLOSURE_REVIEW.md. No production generator or numerical library
is imported. Forty-eight closed rational boxes cover six residues, all h in
[0,1/204], and all true rounding errors. Full-vector replays check the moment
identities and the factor 2 in the quadratic reconstruction independently.
"""
from fractions import Fraction as Q
import argparse
import hashlib
import json
from pathlib import Path
import time


ROOT = Path(__file__).resolve().parents[1]
REVIEW = ROOT / "experiments/v04_all_k/review_2026-10-23"
J = (Q(1999, 8000), Q(2001, 8000))


def require(condition, label):
    if not condition:
        raise ValueError(label)


class Interval:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo, self.hi = Q(lo), Q(lo if hi is None else hi)
        require(self.lo <= self.hi, "reversed interval")

    def __add__(self, other):
        other = interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-interval(other))

    def __rsub__(self, other):
        return interval(other) + (-self)

    def __mul__(self, other):
        other = interval(other)
        products = [a * b for a in (self.lo, self.hi)
                    for b in (other.lo, other.hi)]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = interval(other)
        require(not (other.lo <= 0 <= other.hi), "zero in denominator")
        return self * Interval(1 / other.hi, 1 / other.lo)


def interval(value):
    return value if isinstance(value, Interval) else Interval(value)


def square(value):
    value = interval(value)
    lo2, hi2 = value.lo * value.lo, value.hi * value.hi
    return Interval(0 if value.lo <= 0 <= value.hi else min(lo2, hi2), max(lo2, hi2))


def enclosed_moments(residue, alpha, h, eta):
    """Exact formulas on true grids; interval extension on the parameter box."""
    rho = 1 - alpha
    u0 = Q(1, 2) - (residue % 2) * h / 2
    lam = alpha * (2 * u0 - square(u0))
    u1 = Q(1, 2) + 3 * lam / 4 + eta * h
    nu = alpha * (2 * u1 - square(u1) - lam * (Q(4, 3) - 3 * u1 / 2))
    mu, beta = lam + nu, rho * lam + nu
    delta, theta = Q(residue, 6), Q(6 - residue, 3)
    b = Q(5, 6) + delta * h
    L = b - h
    c = 2 * b - square(b)
    X = 2 - 2 * b + h + mu * theta
    q = theta * h
    gamma = (c - h * X) / (1 - q)
    t1, t2 = Q(residue % 3, 3), Q((2 * residue) % 3, 3)
    v1, v2 = t1 * (1 - t1), t2 * (1 - t2)
    Cg = Q(3, 2) - 2 * theta * h - h * (v1 + v2 / 4)
    segments = [(interval(0), interval(Q(1, 3)), square(rho), beta, Q(1, 2)),
                (interval(Q(1, 3)), u0, square(rho), beta, Q(3, 2)),
                (u0, u1, rho, beta, Q(3, 2)),
                (u1, interval(Q(2, 3)), interval(1), mu, Q(3, 2)),
                (interval(Q(2, 3)), L, interval(1), mu, Q(2))]
    E, Z = interval(0), interval(0)
    for left, right, a, coefficient, g in segments:
        width = right - left
        require(width.lo > 0, "ordered density segments")
        mean = a * (2 - left - right) + (coefficient - gamma) * g
        Z += width * g * mean
        E += width * (square(mean) + square(a) * square(width) / 3)
    # Covariance between the kernel step and the linear / step residual.
    Z += square(h) * (v1 * square(rho) + v2 / 2) - h * (
        v1 * (beta - gamma) + v2 * (mu - gamma) / 4)
    linear_variance = square(h) * (
        square(square(rho)) * u0 + square(rho) * (u1 - u0) + L - u1) / 3
    jump_variance = h * (v1 * square(beta - gamma) + v2 * square(mu - gamma) / 4)
    jump_cross = -2 * square(h) * (
        v1 * square(rho) * (beta - gamma) + v2 * (mu - gamma) / 2)
    E -= linear_variance + jump_variance + jump_cross
    Ahat = 1 + h * Cg / square(1 - q)
    B = -2 * c * Cg / square(1 - q)
    d = -Z / (1 - q)
    f = c * Z / (1 - q)
    inverse_A = h / Ahat
    S3 = 2 * (E - square(d) * inverse_A) * rho + 2 * f - B * d * inverse_A
    S4 = (E * (square(rho) - 1) + 2 * f * (rho - 1)
          - h * X * (Ahat * X + B + 2 * d)
          - square(B + 2 * d * rho) * inverse_A / 4)
    return {"S3": S3, "S4": S4, "E": E, "Z": Z, "Cg": Cg,
            "gamma": gamma, "X": X, "Ahat": Ahat, "B": B,
            "d": d, "f": f, "u0": u0, "u1": u1}


def certify_boxes():
    records = []
    for residue in range(6):
        for j in range(4):
            for k in range(2):
                h = Interval(Q(j, 816), Q(j + 1, 816))
                eta = Interval(Q(-1, 2) + Q(k, 2), Q(-1, 2) + Q(k + 1, 2))
                out = enclosed_moments(residue, Interval(*J), h, eta)
                # The jump cells are wholly in the indicated cutoff regions.
                require(out["u0"].lo > Q(1, 3) + h.hi, "first jump before u0")
                require(out["u1"].hi < Q(2, 3) - h.hi, "second jump after u1")
                require(out["S3"].hi < -Q(3, 25), "uniform derivative margin")
                require(out["S4"].lo > Q(1, 100), "uniform tail value margin")
                records.append({"residue": residue, "h": [str(h.lo), str(h.hi)],
                                "eta": [str(eta.lo), str(eta.hi)],
                                "Hprime_rho_upper": str(out["S3"].hi),
                                "H_rho_lower": str(out["S4"].lo)})
    require(len(records) == 48, "parameter box coverage")
    return records


def cdf(y):
    if y <= Q(1, 3):
        return y / 2
    if y <= Q(2, 3):
        return 3 * y / 2 - Q(1, 3)
    if y <= Q(5, 6):
        return 2 * y - Q(2, 3)
    return Q(1)


def vector_step(p, weights, alpha):
    K, acc, discrepancies = len(p), Q(0), []
    for i, mass in enumerate(p, 1):
        acc += mass
        discrepancies.append(acc - Q(i, K))
    best = max(discrepancies)
    require(best > Q(1, 10**12), "crosscheck does not stop")
    selected = 1 + discrepancies.index(best)
    removed = alpha * sum(p[:selected])
    return [((1 - alpha) * mass if i < selected else mass) + removed * g
            for i, (mass, g) in enumerate(zip(p, weights))], selected


def check_vector_case(K, alpha):
    h, rho, m = Q(1, K), 1 - alpha, (5 * K + 5) // 6
    weights = [cdf(Q(i, K)) - cdf(Q(i - 1, K)) for i in range(1, K + 1)]
    p0 = [Q(2 * K - 2 * i + 1, K * K) for i in range(1, K + 1)]
    p1, s0 = vector_step(p0, weights, alpha)
    p2, s1 = vector_step(p1, weights, alpha)
    require(s0 == K // 2 and p2[m:] == p0[m:], "true prefix and frozen suffix")
    lam = alpha * sum(p0[:s0])
    eta = s1 - K * (Q(1, 2) + 3 * lam / 4)
    require(-Q(1, 2) <= eta <= Q(1, 2), "true rounding error")
    x, q, b = p2[m - 1], weights[m - 1], Q(m, K)
    c = b * (2 - b)
    gamma = (c - x) / (1 - q)
    e = [p2[i] - gamma * weights[i] for i in range(m - 1)]
    require(sum(e) == 0, "E1 vanishes exactly")
    Cg = K * sum(g * g for g in weights[:m - 1])
    E = K * sum(value * value for value in e)
    Z = K * sum(g * value for g, value in zip(weights, e))
    A = K + Cg / (1 - q)**2
    B, d, f = -2 * c * Cg / (1 - q)**2, -Z / (1 - q), c * Z / (1 - q)
    suffix = K * sum((mass - h)**2 for mass in p0[m:])
    C = c * c * Cg / (1 - q)**2 - 2 * c + b + suffix
    U2 = K * sum((mass - h)**2 for mass in p2)
    require(U2 == A*x*x + B*x + C + 2*(d*x + f) + E, "quadratic reconstruction at t=2")
    H2, H1, H0 = E - d*d/A, 2*f - B*d/A, C - U2 - B*B/(4*A)
    require(H2 == (K*(1-q)**2*E + Cg*E - Z*Z)/(A*(1-q)**2),
            "scaled Gram identity on full-vector moments")
    require(H2 >= 0, "Gram sign on crosscheck")
    exact = {"E": E, "Z": Z, "Cg": Cg, "gamma": gamma, "X": x/h,
             "Ahat": h*A, "B": B, "d": d, "f": f,
             "S3": 2*H2*rho + H1, "S4": H2*rho*rho + H1*rho + H0}
    out = enclosed_moments(K % 6, Interval(alpha), Interval(h), Interval(eta))
    for name, expected in exact.items():
        require(out[name].lo == expected == out[name].hi, f"K={K}: exact identity {name}")
    p = p2
    for t in (3, 4):
        p, selected = vector_step(p, weights, alpha)
        require(selected in (m - 1, m), "crosscheck boundary residence")
        z, xt = rho**(t - 2), p[m - 1]
        actual = K * sum((mass - h)**2 for mass in p)
        rebuilt = A*xt*xt + B*xt + C + 2*z*(d*xt + f) + E*z*z
        require(actual == rebuilt, "quadratic reconstruction after a real update")
    return {"K": K, "alpha": str(alpha), "word_head": [s0, s1], "eta": str(eta)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path,
                        help="write a new audit report; default only prints results")
    args = parser.parse_args()
    start = time.monotonic()
    boxes = certify_boxes()
    KS = list(range(33, 75)) + [96, 160, 203, 204, 205, 208, 224, 255, 256,
                              288, 300, 512, 999, 1000, 1001, 2048, 4096]
    checks = [check_vector_case(K, alpha) for K in KS for alpha in (J[0], Q(1, 4), J[1])]
    # Bind the already accepted finite complement; proof-producing modules
    # and mutable checkpoint loaders are not dependencies of that acceptance.
    small_file = REVIEW / "finite_common_audit.json"
    small = json.loads(small_file.read_text(encoding="utf-8"))
    require(small["status"] == "VERIFIED" and not small["failures"], "finite audit status")
    require(set(range(6, 204)) <= set(small["verified_K"]), "finite complement has no holes")
    require(small["common_interval"] == [str(a) for a in J], "finite complement uses J")
    sources = [Path(__file__).resolve(), ROOT / "docs/v04_research/UNIFORM_TAIL_CLOSURE_REVIEW.md",
               ROOT / "docs/v04_research/W0_UNIFORM_BOUND_REVIEW.md",
               ROOT / "docs/v04_research/LARGE_K_PREFIX_LEMMA_REVIEW.md",
               ROOT / "docs/v04_research/D2_ENTRY_REVIEW.md", small_file]
    recorded_hashes = {name.replace("\\", "/"): digest
                       for name, digest in small["source_sha256"].items()}
    for name in ("tools/audit_v04_interval_v3.py", "experiments/v04_w_interval/interval_certs_v3.jsonl",
                 "docs/XO_ALPHA_EXTENSION_PROOF.md", "docs/S3_XO_BOUNDARY_THEOREM_PROOF.md"):
        source = ROOT / name
        require(hashlib.sha256(source.read_bytes()).hexdigest() ==
                recorded_hashes[name], "finite dependency hash: " + name)
        sources.append(source)
    report = {"status": "VERIFIED_UNIFORM_TAIL", "J": [str(a) for a in J],
              "K_min": 204, "time_domain": "all integer t>=3",
              "strict_gap_Ut_minus_U2": "1/100", "Hprime_rho_upper": "-3/25",
              "box_count": len(boxes), "boxes": boxes,
              "full_vector_case_count": len(checks), "full_vector_cases": checks,
              "finite_complement": {"K": [6, 203], "scope": "all alpha in J, all t>=0",
                                    "auditor": "tools/audit_v04_interval_v3.py"},
              "analytic_dependencies": ["true selector prefix lemma on J", "D2 entry and theorem A residence",
                                        "exact bin-variance and covariance identities", "E1=0 by prefix mass",
                                        "Gram sign and completed square", "all-K W0 for t=0,1",
                                        "accepted all-time finite complement"],
              "source_sha256": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sources},
              "elapsed_seconds": round(time.monotonic() - start, 3)}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"VERIFIED: {len(boxes)} exact closed boxes cover all K>=204, alpha in J;")
    print(f"H'(rho)<-3/25, H(rho)>1/100; {len(checks)} exact full-vector checks; "
          f"{report['elapsed_seconds']:.2f}s")


if __name__ == "__main__":
    main()
