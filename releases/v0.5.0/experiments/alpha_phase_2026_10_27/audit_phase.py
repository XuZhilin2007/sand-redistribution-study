"""Independent algebraic/bin-sum checks and analytic prefix constants.

The formal symbols stay formal; point replays are an additional check only.
No experimental status field is treated as a theorem.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
from audit_v040_review import Polynomial as Poly, poly, kernel_integers, integer_step, loss as int_loss
from phase_probe import model, replay
from uniform_probe import enclose, I, ALPHA, BRANCH_ALPHA
from verify_phase import polynomial_states, tail_polynomials, P


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def formal_sums(r, branch):
    h, alpha, eta, gamma = [Poly.variable(i) for i in range(4)]
    rho = 1 - alpha
    u0 = Q(1, 2) - (r % 2) * h / 2
    lam = alpha * (2 * u0 - u0 ** 2)
    u1 = Q(1, 2) + (Q(3, 4) if branch == 'M' else Q(1)) * lam + eta * h
    nu = alpha * (2 * u1 - u1 ** 2 - lam *
                  (Q(4, 3) - 3 * u1 / 2 if branch == 'M' else Q(5, 3) - 2 * u1))
    beta, mu = rho * lam + nu, lam + nu
    b = Q(5, 6) + Q(r, 6) * h
    ell = b - h
    t1, t2, t3 = Q(r % 3, 3), Q(2 * r % 3, 3), Q(5 * r % 6, 6)
    n1, n2, n3 = Q(1, 3) - t1 * h, Q(2, 3) - t2 * h, Q(5, 6) - t3 * h
    first_runs = [
        (h, n1, rho ** 2, beta, Q(1, 2)),
        (n1 + h, n1 + h, rho ** 2, beta, Q(3, 2) - t1),
        (n1 + 2*h, u0, rho ** 2, beta, Q(3, 2))]
    if branch == 'M':
        runs = first_runs + [
            (u0+h, u1, rho, beta, Q(3, 2)),
            (u1+h, n2, poly(1), mu, Q(3, 2)),
            (n2+h, n2+h, poly(1), mu, Q(2)-t2/2),
            (n2+2*h, ell, poly(1), mu, Q(2))]
        segments = [(poly(0), poly(Q(1, 3)), rho**2, beta, Q(1, 2)),
                    (poly(Q(1, 3)), u0, rho**2, beta, Q(3, 2)),
                    (u0, u1, rho, beta, Q(3, 2)),
                    (u1, poly(Q(2, 3)), poly(1), mu, Q(3, 2)),
                    (poly(Q(2, 3)), ell, poly(1), mu, Q(2))]
        a2, coefficient2 = poly(1), mu
    else:
        runs = first_runs + [
            (u0+h, n2, rho, beta, Q(3, 2)),
            (n2+h, n2+h, rho, beta, Q(2)-t2/2),
            (n2+2*h, u1, rho, beta, Q(2)),
            (u1+h, ell, poly(1), mu, Q(2))]
        segments = [(poly(0), poly(Q(1, 3)), rho**2, beta, Q(1, 2)),
                    (poly(Q(1, 3)), u0, rho**2, beta, Q(3, 2)),
                    (u0, poly(Q(2, 3)), rho, beta, Q(3, 2)),
                    (poly(Q(2, 3)), u1, rho, beta, Q(2)),
                    (u1, ell, poly(1), mu, Q(2))]
        a2, coefficient2 = rho, beta

    def sum_runs(rows):
        cg, gp, pp = poly(0), poly(0), poly(0)
        for left, right, a, coefficient, kernel in rows:
            width = right - left + h
            mean = 2 + h - left - right
            sw = width * mean
            sw2 = width * (mean**2 + (width**2 - h**2) / 3)
            cg += width * kernel**2
            gp += kernel * (a * sw + coefficient * kernel * width)
            pp += a**2 * sw2 + 2*a*coefficient*kernel*sw + coefficient**2*kernel**2*width
        return cg, gp, pp

    def integrate(rows, residual=False):
        gp, pp = poly(0), poly(0)
        for left, right, a, coefficient, kernel in rows:
            coefficient = coefficient - gamma if residual else coefficient
            intercept, slope = 2*a + coefficient*kernel, -2*a
            gp += kernel * (intercept*(right-left) + slope*(right**2-left**2)/2)
            pp += (intercept**2*(right-left) + intercept*slope*(right**2-left**2)
                   + slope**2*(right**3-left**3)/3)
        return gp, pp

    cg, gp, pp = sum_runs(runs)
    iz, ie = integrate(segments, residual=True)
    ip, ipp = integrate(segments)
    v1, v2, v3 = t1*(1-t1), t2*(1-t2), t3*(1-t3)
    theta = Q(6-r, 3)
    cg2 = Q(3, 2) - 2*theta*h - h*(v1+v2/4)
    c2 = coefficient2 - gamma
    z2 = iz + h**2*(v1*rho**2+v2*a2/2) - h*(v1*(beta-gamma)+v2*c2/4)
    vlin = h**2*(rho**4*u0+rho**2*(u1-u0)+ell-u1)/3
    vjump = h*(v1*(beta-gamma)**2+v2*c2**2/4)
    vcross = -2*h**2*(v1*rho**2*(beta-gamma)+v2*a2*c2/2)
    e2 = ie - vlin - vjump - vcross
    pp2 = ipp - vlin - h*(v1*beta**2+v2*coefficient2**2/4) + 2*h**2*(v1*rho**2*beta+v2*a2*coefficient2/2)

    # Independent exact sum of p1 across ALL cells, including the support edge.
    one_runs = [(h, n1, rho, lam, Q(1, 2)),
                (n1+h, n1+h, rho, lam, Q(3, 2)-t1),
                (n1+2*h, u0, rho, lam, Q(3, 2)),
                (u0+h, n2, poly(1), lam, Q(3, 2)),
                (n2+h, n2+h, poly(1), lam, Q(2)-t2/2),
                (n2+2*h, n3, poly(1), lam, Q(2)),
                (n3+h, n3+h, poly(1), lam, 2*t3),
                (n3+2*h, poly(1), poly(1), poly(0), Q(0))]
    _, _, pp1_direct = sum_runs(one_runs)
    one_segments = [(poly(0), poly(Q(1, 3)), rho, lam, Q(1, 2)),
                    (poly(Q(1, 3)), u0, rho, lam, Q(3, 2)),
                    (u0, poly(Q(2, 3)), poly(1), lam, Q(3, 2)),
                    (poly(Q(2, 3)), poly(Q(5, 6)), poly(1), lam, Q(2)),
                    (poly(Q(5, 6)), poly(1), poly(1), poly(0), Q(0))]
    _, pp1 = integrate(one_segments)
    pp1 -= h**2*(rho**2*u0+1-u0)/3 + h*lam**2*(v1+v2/4+4*v3)
    pp1 += 2*h**2*lam*(v1*rho+v2/2-2*v3)

    identities = {'Cg': (cg, cg2), 'Z': (gp-gamma*cg, z2),
                  'E': (pp-2*gamma*gp+gamma**2*cg, e2),
                  'prefix_p2_square_norm': (pp, pp2), 'p1_square_norm': (pp1_direct, pp1)}
    result = []
    for name, (left, right) in identities.items():
        require(not (left-right).terms, f'formal identity {branch},{r},{name}')
        result.append({'branch': branch, 'residue': r, 'identity': name,
                       'monomials': len(left.terms), 'sha256': left.digest()})
    if r % 3:
        require((gp-gamma*cg-(z2-2*h**2*(v1*rho**2+v2*a2/2))).terms,
                'control detects reversed covariance')
        require((pp-2*gamma*gp+gamma**2*cg-(ie-vlin-vjump)).terms,
                'control detects omitted mixed term')
    return result


def prefix_constants():
    h = Q(1, 204)
    p = lambda x: x*(2-x)
    fmin = p(Q(1, 2)-h/2)
    lmin, lmax = ALPHA[0]*fmin, ALPHA[1]*Q(3, 4)
    checks = {
        'scanned_first_slope': Q(1, 3)-Q(4, 3)*ALPHA[1]+lmin/2,
        'scanned_middle_slope': ALPHA[0]*(Q(3, 2)*fmin-1),
        'middle_peak_above_u0': Q(3, 4)*lmin-h/2,
        'middle_peak_below_2thirds': Q(2, 3)-(Q(1, 2)+Q(3, 4)*lmax+h/2),
        'upper_peak_above_2thirds': Q(1, 2)+lmin-h/2-Q(2, 3),
        'upper_peak_below_5sixths': Q(5, 6)-(Q(1, 2)+lmax+h/2),
        'middle_dominates_scanned': (Q(3, 4)*lmin)**2-h*h/4,
        'upper_dominates_suffix': (Q(5, 6)-Q(1, 2)-lmax)**2-h*h/4,
        'M_dominates_when_alpha_at_most_quarter': Q(3, 16)*(4-21*Q(3, 16))/48-h*h/4,
        'U_dominates_when_alpha_at_least_13over50':
            (Q(13, 50)*fmin)*(21*Q(13, 50)*fmin-4)/48-h*h/4,
        'lambda_above_peak_gap_vertex': lmin-Q(2, 21),
    }
    require(all(v > 0 for v in checks.values()), 'uniform prefix constants')
    require(BRANCH_ALPHA['M'][1] >= BRANCH_ALPHA['U'][0], 'two branch domains overlap')
    return {k: str(v) for k, v in checks.items()}


def vector_checks():
    records = []
    for k in list(range(204, 216)) + [224, 288, 300, 999, 1000, 1001, 2048]:
        for alpha in (ALPHA[0], Q(1, 4), Q(16, 63), Q(13, 50), ALPHA[1]):
            ps, word, uu = replay(k, alpha, 2)
            p0, g, h = model(k, True)
            m = (5*k+5)//6
            c, q = Q(m, k)*(2-Q(m, k)), g[m-1]
            u1, lam = Q(word[1], k), alpha*sum(p0[:k//2])
            branch = 'M' if u1 < Q(2, 3) else 'U'
            eta = (u1-Q(1, 2)-(Q(3, 4) if branch == 'M' else 1)*lam)/h
            require(-Q(1, 2) <= eta <= Q(1, 2), 'true prefix enclosure')
            require(BRANCH_ALPHA[branch][0] <= alpha <= BRANCH_ALPHA[branch][1], 'true branch domain')
            values = enclose(k%6, I(alpha), I(h), I(eta), branch)
            gamma = (c-ps[2][m-1])/(1-q)
            error = [ps[2][i]-gamma*g[i] for i in range(m-1)]
            expected = {'U1': uu[1], 'U2': uu[2], 'Cg': k*sum(v*v for v in g[:m-1]),
                        'Z': k*sum(v*w for v, w in zip(g, error)), 'E': k*sum(v*v for v in error)}
            for name, value in expected.items():
                require(values[name].lo == values[name].hi == value, f'full vector identity {name}')
            polys, gg, hh = polynomial_states(k, word)
            h2, h1, h0 = tail_polynomials(k, polys, gg)
            rho = 1-alpha
            require(values['Hrho'].lo == P.value(h2, alpha)*rho*rho+P.value(h1, alpha)*rho+P.value(h0, alpha), 'Hrho identity')
            require(values['-Hprime'].lo == -(2*P.value(h2, alpha)*rho+P.value(h1, alpha)), 'Hprime identity')
            records.append({'K': k, 'alpha': str(alpha), 'branch': branch, 'eta': str(eta), 'word': word})
    # Entirely different integer/common-denominator implementation, with
    # kernel weights formed by intersection lengths rather than CDF values.
    for k in (6, 7, 9, 10, 13, 16, 32, 96, 204, 205):
        for alpha in (Q(6, 25), Q(1, 4), Q(13, 40), Q(9, 10)):
            ps, word, uu = replay(k, alpha, 6)
            nums, den = list(range(2*k-1, 0, -2)), k*k
            weights = kernel_integers(k)
            for t in range(6):
                require(int_loss(nums, den) == uu[t], 'independent integer loss')
                nums, den, s, ties = integer_step(nums, den, weights, alpha)
                require(s == word[t], 'independent integer selector')
                require([Q(n, den) for n in nums] == ps[t+1], 'independent integer state')
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, default=Path(__file__).with_name('algebra_audit.json'))
    args = parser.parse_args()
    started = time.monotonic()
    constants = prefix_constants()
    identities = [record for branch in ('M', 'U') for r in range(6) for record in formal_sums(r, branch)]
    points = vector_checks()
    report = {'status': 'VERIFIED_ALGEBRA_AND_PREFIX', 'prefix_constants': constants,
              'formal_identities': identities, 'formal_identity_count': len(identities),
              'exact_vector_cases': points, 'integer_vector_cases': 40,
              'elapsed_seconds': time.monotonic()-started,
              'source_sha256': {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [Path(__file__), Path(__file__).with_name('uniform_probe.py'),
                                          Path(__file__).with_name('phase_probe.py'), ROOT/'tools/audit_v040_review.py']}}
    args.report.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('PASS', len(identities), 'formal identities;', len(points), 'moment cases; 40 integer replays;',
          report['elapsed_seconds'], 'seconds')


if __name__ == '__main__':
    main()
