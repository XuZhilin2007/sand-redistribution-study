"""Post-release audit of v0.4.0; this is not part of its frozen proof.

Standard library only. Reconstructs scaled tail moments from finite arithmetic
progressions of bin averages, then compares polynomial coefficients with an
independent continuous-integral/variance construction. All four indeterminates
(h, alpha, eta, gamma) remain formal. No Sand modules are imported.

Also replays selected trajectories using common-denominator integer masses,
and checks published counterexamples and the two-slope conjugacy at the tie.
Finite trajectory checks are explicitly not an infinite-time proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from functools import reduce
import hashlib
import json
from math import comb, gcd, lcm
from pathlib import Path
import platform
import time


ROOT = Path(__file__).resolve().parents[1]
J = (Q(1999, 8000), Q(2001, 8000))


def require(condition, message):
    if not condition:
        raise ValueError(message)


class Polynomial:
    """Sparse Q[h,alpha,eta,gamma], with coefficientwise exact equality."""

    def __init__(self, value=0):
        if isinstance(value, dict):
            self.terms = {m: Q(c) for m, c in value.items() if c}
        else:
            self.terms = {(0, 0, 0, 0): Q(value)} if value else {}

    @classmethod
    def variable(cls, i):
        m = [0] * 4
        m[i] = 1
        return cls({tuple(m): Q(1)})

    def __add__(self, other):
        other = poly(other)
        terms = dict(self.terms)
        for m, c in other.terms.items():
            terms[m] = terms.get(m, Q(0)) + c
        return Polynomial(terms)

    __radd__ = __add__

    def __neg__(self):
        return Polynomial({m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-poly(other))

    def __rsub__(self, other):
        return poly(other) + (-self)

    def __mul__(self, other):
        other = poly(other)
        terms = {}
        for m, c in self.terms.items():
            for n, d in other.terms.items():
                key = tuple(x + y for x, y in zip(m, n))
                terms[key] = terms.get(key, Q(0)) + c * d
        return Polynomial(terms)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        return self * (1 / Q(scalar))

    def __pow__(self, exponent):
        require(type(exponent) is int and exponent >= 0, "nonnegative polynomial power")
        result = Polynomial(1)
        for _ in range(exponent):
            result = result * self
        return result

    def evaluate(self, values):
        return sum((c * reduce(lambda x, y: x * y,
                              (v**n for v, n in zip(values, m)), Q(1))
                    for m, c in self.terms.items()), Q(0))

    def digest(self):
        payload = [[list(m), str(c)] for m, c in sorted(self.terms.items())]
        return hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()


def poly(value):
    return value if isinstance(value, Polynomial) else Polynomial(value)


def formal_moments(residue):
    h, alpha, eta, gamma = [Polynomial.variable(i) for i in range(4)]
    rho = 1 - alpha
    u0 = Q(1, 2) - (residue % 2) * h / 2
    lam = alpha * (2 * u0 - u0**2)
    u1 = Q(1, 2) + 3 * lam / 4 + eta * h
    nu = alpha * (2 * u1 - u1**2 - lam * (Q(4, 3) - 3 * u1 / 2))
    beta, mu = rho * lam + nu, lam + nu
    ell = Q(5, 6) + (Q(residue, 6) - 1) * h
    t1, t2 = Q(residue % 3, 3), Q((2 * residue) % 3, 3)
    n1h, n2h = Q(1, 3) - t1 * h, Q(2, 3) - t2 * h

    # Inclusive endpoints a=i_min*h, b=i_max*h. Each cell's initial
    # average density is w_i=2+h-2*i*h. The two designated singleton
    # cells also work when t=0: their kernel is the full right value.
    runs = [
        (h, n1h, rho**2, beta, Q(1, 2)),
        (n1h + h, n1h + h, rho**2, beta, Q(3, 2) - t1),
        (n1h + 2*h, u0, rho**2, beta, Q(3, 2)),
        (u0 + h, u1, rho, beta, Q(3, 2)),
        (u1 + h, n2h, poly(1), mu, Q(3, 2)),
        (n2h + h, n2h + h, poly(1), mu, Q(2) - t2 / 2),
        (n2h + 2*h, ell, poly(1), mu, Q(2)),
    ]
    cg, gp, pp = poly(0), poly(0), poly(0)
    for a, b, A, C, k in runs:
        width = b - a + h  # h times the number of cells
        mean = 2 + h - a - b
        sum_w = width * mean
        sum_w2 = width * (mean**2 + (width**2 - h**2) / 3)
        cg += width * k*k
        gp += k * (A * sum_w + C * k * width)
        pp += A*A * sum_w2 + 2*A*C*k * sum_w + C*C*k*k * width
    direct = {"Cg": cg, "Z": gp - gamma*cg,
              "E": pp - 2*gamma*gp + gamma*gamma*cg}

    # A separate route: integrate an affine residual on each physical
    # segment, then remove the exact within-cell variances/covariances.
    segments = [
        (poly(0), poly(Q(1, 3)), rho**2, beta, Q(1, 2)),
        (poly(Q(1, 3)), u0, rho**2, beta, Q(3, 2)),
        (u0, u1, rho, beta, Q(3, 2)),
        (u1, poly(Q(2, 3)), poly(1), mu, Q(3, 2)),
        (poly(Q(2, 3)), ell, poly(1), mu, Q(2)),
    ]
    ie, iz = poly(0), poly(0)
    for a, b, A, C, k in segments:
        intercept, slope = 2*A + (C-gamma)*k, -2*A
        iz += k * (intercept*(b-a) + slope*(b*b-a*a)/2)
        ie += (intercept*intercept*(b-a)
               + intercept*slope*(b*b-a*a)
               + slope*slope*(b**3-a**3)/3)
    v1, v2 = t1*(1-t1), t2*(1-t2)
    theta = Q(6-residue, 3)
    cg_integral = Q(3, 2) - 2*theta*h - h*(v1+v2/4)
    z_integral = iz + h*h*(v1*rho*rho+v2/2) - h*(v1*(beta-gamma)+v2*(mu-gamma)/4)
    linear = h*h * (rho**4*u0+rho*rho*(u1-u0)+ell-u1)/3
    jump = h * (v1*(beta-gamma)**2+v2*(mu-gamma)**2/4)
    cross = -2*h*h * (v1*rho*rho*(beta-gamma)+v2*(mu-gamma)/2)
    integral = {"Cg": cg_integral, "Z": z_integral, "E": ie-linear-jump-cross}
    records = []
    for name in direct:
        difference = direct[name] - integral[name]
        require(not difference.terms, f"residue {residue}: polynomial identity {name}")
        if name == "Z" and residue % 3:
            wrong_covariance = z_integral - 2*h*h*(v1*rho*rho+v2/2)
            require((direct[name]-wrong_covariance).terms,
                    "covariance sign error must be detected symbolically")
        if name == "E":
            require((direct[name]-(ie-linear-jump)).terms or residue % 3 == 0,
                    "omitted mixed term must be detected symbolically")
        records.append({"residue": residue, "moment": name,
                        "monomial_count": len(direct[name].terms),
                        "polynomial_sha256": direct[name].digest(),
                        "difference_monomials": len(difference.terms)})
    return direct, records


def kernel_integers(k):
    # Compute the kernel by lengths of intersection with density segments,
    # without evaluating the CDF used by the released implementations.
    knots = (Q(0), Q(1, 3), Q(2, 3), Q(5, 6), Q(1))
    densities = (Q(1, 2), Q(3, 2), Q(2), Q(0))
    weights = []
    for i in range(k):
        left, right = Q(i, k), Q(i+1, k)
        mass = sum((max(Q(0), min(right, b)-max(left, a))*d
                    for a, b, d in zip(knots, knots[1:], densities)), Q(0))
        scaled = mass * (6*k)
        require(scaled.denominator == 1, "integer kernel scale")
        weights.append(scaled.numerator)
    require(sum(weights) == 6*k and min(weights) >= 0, "kernel normalization")
    return weights


def integer_step(nums, denominator, weights, alpha):
    k, running, best, selected, at_selected = len(nums), 0, None, 0, 0
    ties = []
    for i, n in enumerate(nums, 1):
        running += n
        score = k*running - i*denominator
        if best is None or score > best:
            best, selected, at_selected, ties = score, i, running, [i]
        elif score == best:
            ties.append(i)
    require(best * 10**12 > k*denominator, "trajectory stopped")
    a, b, scale = alpha.numerator, alpha.denominator, 6*k
    result = [(b-a if i < selected else b)*n*scale + a*at_selected*w
              for i, (n, w) in enumerate(zip(nums, weights))]
    den = b*denominator*scale
    common = reduce(gcd, result, den)
    result, den = [n//common for n in result], den//common
    require(sum(result) == den and min(result) >= 0, "mass/nonnegativity")
    return result, den, selected, ties


def loss(nums, denominator):
    return Q(len(nums)*sum(n*n for n in nums), denominator*denominator)-1


def replay(k, alpha, horizon=60):
    nums = list(range(2*k-1, 0, -2))
    den, weights = k*k, kernel_integers(k)
    energies, word, ties, state2 = [], [], [], None
    m, h = (5*k+5)//6, Q(1, k)
    c = Q(m, k) * (2-Q(m, k))
    q = Q(weights[m-1], 6*k)
    for t in range(horizon+1):
        energies.append(loss(nums, den))
        if t == 2:
            state2 = [Q(n, den) for n in nums]
        if t == horizon:
            break
        previous_x = Q(nums[m-1], den)
        nums, den, selected, maxima = integer_step(nums, den, weights, alpha)
        word.append(selected)
        if len(maxima) > 1:
            ties.append({"t": t, "maxima": maxima})
        require(Q(sum(nums[:m]), den) == c, "conserved prefix")
        if t >= 2:
            require(selected in (m-1, m), "boundary residence in finite replay")
            expected = ((1-alpha)*previous_x+alpha*c*q if previous_x > h
                        else (1-alpha*q)*previous_x+alpha*c*q)
            require(Q(nums[m-1], den) == expected, "scalar recurrence")
    require(min(energies) == energies[2], "finite two-update optimum")
    record = {"K": k, "alpha": str(alpha), "horizon": horizon,
              "word_prefix": word[:6], "ties": ties,
              "minimum_times": [i for i, u in enumerate(energies) if u == energies[2]],
              "U0_minus_U2": str(energies[0]-energies[2]),
              "U1_minus_U2": str(energies[1]-energies[2]),
              "min_Ut_minus_U2_for_3_to_horizon": str(min(energies[3:])-energies[2])}
    if k >= 33:
        p0 = [Q(n, k*k) for n in range(2*k-1, 0, -2)]
        lam = alpha*sum(p0[:word[0]])
        eta = word[1] - k*(Q(1, 2)+3*lam/4)
        require(-Q(1, 2) <= eta <= Q(1, 2), "rounding enclosure")
        cg = k*sum(Q(w, 6*k)**2 for w in weights[:m-1])
        gamma = (c-state2[m-1])/(1-q)
        error = [p-gamma*Q(w, 6*k) for p, w in zip(state2[:m-1], weights)]
        observed = {"Cg": cg,
                    "Z": k*sum(Q(w, 6*k)*e for w, e in zip(weights, error)),
                    "E": k*sum(e*e for e in error)}
        for name, expression in MOMENTS[k % 6].items():
            require(expression.evaluate((h, alpha, eta, gamma)) == observed[name],
                    f"bin-sum symbolic formula vs integer replay: {name}")
        record["symbolic_moments_match_integer_state"] = True
    return record


def conjugacy_checks():
    records = []
    for k in (9, 10, 11, 21, 65, 195, 479):
        for alpha in (Q(1, 8), Q(1, 4), Q(7, 24)):
            h, m = Q(1, k), (5*k+5)//6
            c, q = Q(m, k)*(2-Q(m, k)), Q(6-k%6, 3*k)
            L = c*q
            A = lambda x: (1-alpha)*x+alpha*L
            B = lambda x: (1-alpha*q)*x+alpha*L
            lower, upper = A(h), B(h)
            width = upper-lower
            eta = (upper-h)/width
            require(eta == q*(c-h)/((1-q)*h), "threshold dictionary")
            lam, mu = 1-alpha, (1-alpha*q)/(1-alpha)
            delta = 1-lam*eta
            require(0 < lam < 1 and 1-lam < delta < 1 and lam*mu < 1,
                    "Laurent-Nogueira hypotheses")
            for z in (Q(0), eta/2, eta, (1+eta)/2, Q(999, 1000)):
                x = upper-width*z
                image = (upper-(A(x) if x > h else B(x)))/width
                standard = lam*z+delta if z < eta else mu*(lam*z+delta-1)
                require(image == standard, "conjugacy including exact tie")
            require((upper-B(h))/width == 0, "right branch at threshold")
            records.append({"K": k, "alpha": str(alpha), "eta": str(eta),
                            "lambda": str(lam), "mu": str(mu), "delta": str(delta)})
    return records


def counterexamples():
    records = []
    cases = [(9, Q(1, 5), 3, Q(-83303117, 5978711250)),
             (6, Q(7, 24), 1, Q(-161, 11943936))]
    for k, alpha, t, expected in cases:
        nums, den, weights = list(range(2*k-1, 0, -2)), k*k, kernel_integers(k)
        energies, word = [], []
        for step in range(4):
            energies.append(loss(nums, den))
            if step < 3:
                nums, den, selected, _ = integer_step(nums, den, weights, alpha)
                word.append(selected)
        difference = energies[t]-energies[2]
        require(difference == expected < 0, "published exact counterexample")
        records.append({"K": k, "alpha": str(alpha), "t": t,
                        "Ut_minus_U2": str(difference), "word": word})
    return records


def historical_m0_checks():
    # General fixed-zone M0, not adaptive XO. This legal near-biased q
    # has a continuous interior vertex, but a coarse time step skips it.
    q = [Q(17, 100)]*3 + [Q(49, 300)]*3
    near, far = sum(q[:3]), sum(q[3:])
    sn, sf = sum(p*p for p in q[:3]), sum(p*p for p in q[3:])
    curvature = sn + sf*near*near/(far*far)
    vertex = near*sf/(far*far*sn+near*near*sf)
    require(near+far == 1 and sn/near > sf/far and vertex == Q(50, 51),
            "M0 counterexample satisfies historical hypotheses")

    def m0loss(z):
        return sn*z*z + sf*((1-near*z)/far)**2 - Q(1, 6)

    alpha = Q(1, 4)
    rho = 1-alpha*far
    direct = [(1-alpha)*p+alpha*near*p if i < 3 else p+alpha*near*p
              for i, p in enumerate(q)]
    observed = sum((p-Q(1, 6))**2 for p in direct)
    require(observed == m0loss(rho), "M0 direct first update")
    require(0 < rho < vertex < 1 and m0loss(rho) > m0loss(1),
            "all later z lie beyond the vertex, so t=0 uniquely minimizes")
    small_rho = 1-Q(1, 100)*far
    require(m0loss(small_rho) < m0loss(1), "discrete improvement depends on alpha")
    # The exact comparison factors as A*(rho-1)*(rho+1-2*z_star).
    require(m0loss(rho)-m0loss(1) == curvature*(rho-1)*(rho+1-2*vertex),
            "corrected M0 discrete criterion")

    # Exact enumeration of two independent draws, rather than a simulated
    # noise estimate, distinguishes the two claimed expectation formulas.
    n, p = 2, Q(3, 4)
    expected = sum((comb(n, j)*p**j*(1-p)**(n-j)
                    * 2*(Q(j, n)-Q(1, 2))**2 for j in range(n+1)), Q(0))
    udet = 2*(p-Q(1, 2))**2
    corrected = udet+(1-p*p-(1-p)**2)/n
    historical = udet+Q(1, 2)/n
    require(expected == corrected == Q(5, 16) and expected != historical,
            "exact multinomial expectation correction")
    return {"fixed_zone_counterexample": {
                "K": 6, "q": [str(p) for p in q], "near_bins": 3,
                "Q": str(near), "cn": str(sn/near), "cf": str(sf/far),
                "z_star": str(vertex), "alpha": str(alpha), "rho": str(rho),
                "U0": str(m0loss(1)), "U1": str(m0loss(rho)),
                "U1_minus_U0": str(m0loss(rho)-m0loss(1)),
                "conclusion": "t=0 is the unique optimum over all integer t>=0",
                "alpha_1_over_100_U1_minus_U0": str(m0loss(small_rho)-m0loss(1)),
                "correct_strict_improvement_criterion": "z_star < (1+rho)/2"},
            "multinomial_expectation": {
                "N": n, "p": ["3/4", "1/4"], "Udet": str(udet),
                "exact_enumeration": str(expected), "historical_formula": str(historical),
                "correct_formula": "E[Uhat] = Udet + (1-sum(p_i^2))/N"}}


def scope_counterexamples():
    records = []
    k, alpha = 9, Q(1, 4)
    p0 = [Q(n, k*k) for n in range(2*k-1, 0, -2)]
    for perturb in (False, True):
        p = list(p0)
        m = (5*k+5)//6
        if perturb:
            p[m-4] -= Q(1, 10)
            p[m-2] += Q(1, 10)
        require(min(p) >= 0 and sum(p) == 1 and sum(p[:m]) == sum(p0[:m]),
                "legal prefix-mass-preserving perturbation")
        den = lcm(*(a.denominator for a in p))
        nums, weights = [int(a*den) for a in p], kernel_integers(k)
        us, tvs = [], []
        for t in range(4):
            us.append(loss(nums, den))
            tvs.append(sum(abs(Q(a, den)-Q(1, k)) for a in nums)/2)
            if t < 3:
                nums, den, _, _ = integer_step(nums, den, weights, alpha)
        if perturb:
            require(us[3] < us[2], "noncanonical state destroys two-update optimum")
            records.append({"case": "K9 P3 delta=1/10", "L1_change": "1/5",
                            "U3_minus_U2": str(us[3]-us[2])})
        else:
            require(tvs[3] < tvs[2], "total variation differs from quadratic metric")
            records.append({"case": "K9 canonical total variation",
                            "TV3_minus_TV2": str(tvs[3]-tvs[2])})
    return records


MOMENTS = {}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    identities = []
    # A low-degree independent check of the sparse arithmetic, with exact
    # evaluation of multiplication/power and distributivity at rational points.
    x, y, z, w = [Polynomial.variable(i) for i in range(4)]
    sample = (Q(1, 17), Q(1, 4), Q(-1, 3), Q(7, 5))
    p = (x+2*y-z*w)**3
    require(p.evaluate(sample) == (sample[0]+2*sample[1]-sample[2]*sample[3])**3,
            "sparse arithmetic control")
    for residue in range(6):
        MOMENTS[residue], records = formal_moments(residue)
        identities.extend(records)
    print("PASS: 18 formal moment identities in Q[h,alpha,eta,gamma]", flush=True)
    ks = [*range(6, 19), 32, 33, 34, 35, 36, 37, 38, 96, 160, 203, 204, 205, 224, 288]
    cases = [(k, alpha) for k in ks for alpha in (J[0], Q(1, 4), J[1])]
    cases += [(k, Q(1, 4)+epsilon) for k in (16, 32, 96, 160, 224, 288)
              for epsilon in (Q(-1, 10**12), Q(1, 10**12))]
    trajectories = [replay(k, alpha) for k, alpha in cases]
    print(f"PASS: {len(trajectories)} integer-mass trajectories through t=60", flush=True)
    conjugacy = conjugacy_checks()
    counters = counterexamples()
    m0 = historical_m0_checks()
    scope = scope_counterexamples()
    source = Path(__file__).resolve()
    dependencies = [source, ROOT / "docs/V040_REGISTRATION.md",
                    ROOT / "docs/v04_research/UNIFORM_TAIL_CLOSURE_REVIEW.md",
                    ROOT / "docs/XO_ALPHA_EXTENSION_PROOF.md",
                    ROOT / "docs/M0_GENERAL_Q.md", ROOT / "docs/M0_MECHANISM.md"]
    report = {"review_date": "2026-10-27", "status": "PASS",
              "python": platform.python_version(),
              "symbolic_identity_count": len(identities), "symbolic_identities": identities,
              "formal_variables": ["h", "alpha", "eta", "gamma"],
              "trajectory_case_count": len(trajectories), "trajectories": trajectories,
              "conjugacy_case_count": len(conjugacy), "conjugacy": conjugacy,
              "counterexamples": counters,
              "historical_m0_corrections": m0, "scope_counterexamples": scope,
              "limitations": ["Finite trajectories do not prove all-time optimality.",
                              "Polynomial identities verify the stated moment formulas, not their sign on the parameter domain.",
                              "Boundary residence and real-grid alignment remain analytic proof obligations.",
                              "This Python implementation is not a proof-assistant formalization."],
              "source_sha256": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in dependencies},
              "elapsed_seconds": round(time.monotonic()-started, 3)}
    if args.report:
        path = args.report.resolve()
        require(ROOT in path.parents and path != source, "report must be a separate workspace file")
        require(not path.exists(), "refuse to overwrite an existing evidence record")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(f"PASS: {len(conjugacy)} conjugacy cases; 2 exact counterexamples; "
          f"{report['elapsed_seconds']:.2f}s", flush=True)
    print("CONFIRMED: historical M0 discrete-optimum and sampling-expectation corrections", flush=True)


if __name__ == "__main__":
    main()
