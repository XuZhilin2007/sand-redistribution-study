"""Uniform extension: two possible first-step peaks, exact boxes.

Analytic premises, exact bin-sum audit, and finite complement are separate
obligations. Uses the frozen Fraction interval arithmetic as a library.
"""
from __future__ import annotations

from fractions import Fraction as Q
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
from audit_v04_uniform_tail import Interval as I, square as sq

ALPHA = (Q(123, 500), Q(7, 24))
BRANCH_ALPHA = {'M': (ALPHA[0], Q(13, 50)), 'U': (Q(1, 4), ALPHA[1])}
FLOORS = {'U0-U2': Q(3, 20), 'U1-U2': Q(1, 1000), '-Hprime': Q(1, 10),
          'Hrho': Q(1, 200), 'entry_margin': Q(3, 100),
          'upper_slope': Q(0), 'middle_slope': Q(0), 'transition_slope': Q(0)}


def mean_square(a, coefficient, kernel, left, right):
    width = right - left
    mean = a * (2 - left - right) + coefficient * kernel
    return width * (sq(mean) + sq(a) * sq(width) / 3)


def enclose(r, alpha, h, eta, branch):
    rho = 1 - alpha
    u0 = Q(1, 2) - (r % 2) * h / 2
    f0 = 2 * u0 - sq(u0)
    lam = alpha * f0
    u1 = Q(1, 2) + (Q(3, 4) if branch == 'M' else Q(1)) * lam + eta * h
    complement_G = Q(4, 3) - 3 * u1 / 2 if branch == 'M' else Q(5, 3) - 2 * u1
    nu = alpha * (2 * u1 - sq(u1) - lam * complement_G)
    mu, beta = lam + nu, rho * lam + nu
    theta, delta = Q(6 - r, 3), Q(r, 6)
    b = Q(5, 6) + delta * h
    L = b - h
    c, q = 2 * b - sq(b), theta * h
    X = 2 - 2 * b + h + mu * theta
    gamma = (c - h * X) / (1 - q)
    t1, t2, t3 = Q(r % 3, 3), Q(2 * r % 3, 3), Q(5 * r % 6, 6)
    v1, v2, v3 = t1 * (1 - t1), t2 * (1 - t2), t3 * (1 - t3)
    Cg = Q(3, 2) - 2 * theta * h - h * (v1 + v2 / 4)
    if branch == 'M':
        segments = [(I(0), I(Q(1, 3)), sq(rho), beta, Q(1, 2)),
                    (I(Q(1, 3)), u0, sq(rho), beta, Q(3, 2)),
                    (u0, u1, rho, beta, Q(3, 2)),
                    (u1, I(Q(2, 3)), I(1), mu, Q(3, 2)),
                    (I(Q(2, 3)), L, I(1), mu, Q(2))]
        a2, c2 = I(1), mu - gamma
    else:
        segments = [(I(0), I(Q(1, 3)), sq(rho), beta, Q(1, 2)),
                    (I(Q(1, 3)), u0, sq(rho), beta, Q(3, 2)),
                    (u0, I(Q(2, 3)), rho, beta, Q(3, 2)),
                    (I(Q(2, 3)), u1, rho, beta, Q(2)),
                    (u1, L, I(1), mu, Q(2))]
        a2, c2 = rho, beta - gamma
    Z, E = I(0), I(0)
    for left, right, a, coefficient, kernel in segments:
        if (right - left).lo <= 0:
            raise ValueError('unordered density segments')
        width = right - left
        mean = a * (2 - left - right) + (coefficient - gamma) * kernel
        Z += width * kernel * mean
        E += mean_square(a, coefficient - gamma, kernel, left, right)
    Z += sq(h) * (v1 * sq(rho) + v2 * a2 / 2) - h * (v1 * (beta - gamma) + v2 * c2 / 4)
    Vlin = sq(h) * (sq(sq(rho)) * u0 + sq(rho) * (u1 - u0) + L - u1) / 3
    Vjump = h * (v1 * sq(beta - gamma) + v2 * sq(c2) / 4)
    Vcross = -2 * sq(h) * (v1 * sq(rho) * (beta - gamma) + v2 * a2 * c2 / 2)
    E -= Vlin + Vjump + Vcross
    Ahat = 1 + h * Cg / sq(1 - q)
    B = -2 * c * Cg / sq(1 - q)
    d, f = -Z / (1 - q), c * Z / (1 - q)
    inv_A = h / Ahat
    hp = 2 * (E - sq(d) * inv_A) * rho + 2 * f - B * d * inv_A
    hw = E * (sq(rho) - 1) + 2 * f * (rho - 1) - h * X * (Ahat * X + B + 2 * d) - sq(B + 2 * d * rho) * inv_A / 4

    # Exact U2 from the complete vector norm, including the frozen suffix.
    suffix_cubes = (1 - b) * sq(1 - b)
    suffix_norm = Q(4, 3) * suffix_cubes - (1 - b) * sq(h) / 3
    # Integrate f2 directly for W0. Reusing gamma,E,Z here is exact but
    # loses too much dependency information under natural interval extension.
    U2 = sum((mean_square(a, coefficient, kernel, left, right)
              for left, right, a, coefficient, kernel in segments), I(0))
    coeff2 = mu if branch == 'M' else beta
    U2 -= Vlin + h * (v1 * sq(beta) + v2 * sq(coeff2) / 4)
    U2 += 2 * sq(h) * (v1 * sq(rho) * beta + v2 * a2 * coeff2 / 2)
    U2 += h * sq(X) + suffix_norm - 1
    U0 = (1 - sq(h)) / 3
    # Exact U1, including all three mixed kernel cells and their covariance.
    one_segments = [(I(0), I(Q(1, 3)), rho, lam, Q(1, 2)),
                    (I(Q(1, 3)), u0, rho, lam, Q(3, 2)),
                    (u0, I(Q(2, 3)), I(1), lam, Q(3, 2)),
                    (I(Q(2, 3)), I(Q(5, 6)), I(1), lam, Q(2)),
                    (I(Q(5, 6)), I(1), I(1), I(0), Q(0))]
    U1 = sum((mean_square(a, coeff, kern, left, right)
              for left, right, a, coeff, kern in one_segments), I(0)) - 1
    U1 -= sq(h) * (sq(rho) * u0 + 1 - u0) / 3
    U1 -= h * sq(lam) * (v1 + v2 / 4 + 4 * v3)
    U1 += 2 * sq(h) * lam * (v1 * rho + v2 / 2 - 2 * v3)

    # D2 entry: all lower two-thirds coordinates below a boundary value;
    # both pieces of the upper third are strictly increasing.
    a_lo, a_hi = rho.lo, rho.hi
    beta_lo, beta_hi = beta.lo, beta.hi
    Amax = (2 * a_hi ** 2 + beta_hi / 2 - 1) ** 2 / (4 * a_lo ** 2)
    Bmax = -beta_lo / 3 + (2 * a_hi ** 2 + 3 * beta_hi / 2 - 1) ** 2 / (4 * a_lo ** 2)
    Cmax = -4 * a_lo * lam.lo / 3 - nu.lo / 3 + (2 * a_hi + 3 * beta_hi / 2 - 1) ** 2 / (4 * a_lo)
    Emax = Q(2, 9) - mu.lo / 3
    boundary_lower = Q(5, 36) - (2 * mu.hi - Q(1, 3)) * h.hi
    entry_margin = boundary_lower - max(Amax, Bmax, Cmax, Emax)
    upper_slope = 2 * mu.lo - Q(2, 3)
    middle_slope = 1 + Q(3, 2) * mu.lo - Q(4, 3)
    transition_slope = (2 * rho * (1 - u1) + 2 * beta - 1).lo if branch == 'U' else Q(1)
    return {'U0-U2': U0 - U2, 'U1-U2': U1 - U2, '-Hprime': -hp, 'Hrho': hw,
            'entry_margin': I(entry_margin), 'upper_slope': I(upper_slope),
            'middle_slope': I(middle_slope), 'transition_slope': I(transition_slope),
            'Cg': Cg, 'Z': Z, 'E': E, 'U1': U1, 'U2': U2, 'u0': u0, 'u1': u1}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alpha-parts', type=int, default=64)
    parser.add_argument('--h-parts', type=int, default=8)
    parser.add_argument('--K-min', type=int, default=204)
    parser.add_argument('--report', type=Path, default=Path(__file__).with_name('uniform_boxes.json'))
    args = parser.parse_args()
    started = time.monotonic()
    minima = {}
    failures = []
    count = 0
    leaves = []
    def visit(branch, r, aa, hh, eta, depth=0):
        nonlocal count
        count += 1
        result = enclose(r, aa, hh, eta, branch)
        names = ('U0-U2', 'U1-U2', '-Hprime', 'Hrho', 'entry_margin',
                 'upper_slope', 'middle_slope', 'transition_slope')
        failed = [name for name in names if result[name].lo <= FLOORS[name]]
        key = [branch, r, [str(aa.lo), str(aa.hi)], [str(hh.lo), str(hh.hi)],
               [str(eta.lo), str(eta.hi)]]
        if failed and depth < 12:
            intervals = [aa, hh, eta]
            scores = [aa.hi - aa.lo, hh.hi - hh.lo, (eta.hi - eta.lo) * hh.hi]
            dim = max(range(3), key=lambda j: scores[j])
            mid = (intervals[dim].lo + intervals[dim].hi) / 2
            for part in (I(intervals[dim].lo, mid), I(mid, intervals[dim].hi)):
                child = list(intervals)
                child[dim] = part
                visit(branch, r, *child, depth + 1)
            return
        for name in names:
            value = result[name].lo
            if name not in minima or value < Q(minima[name]['lower']):
                minima[name] = {'lower': str(value), 'approximation': float(value), 'box': key}
        if failed:
            failures.append({'obligations': failed, 'box': key})
        leaves.append(key)
    for branch in ('M', 'U'):
        alpha_domain = BRANCH_ALPHA[branch]
        for r in range(6):
            for ai in range(args.alpha_parts):
                aa = I(alpha_domain[0] + (alpha_domain[1] - alpha_domain[0]) * Q(ai, args.alpha_parts),
                       alpha_domain[0] + (alpha_domain[1] - alpha_domain[0]) * Q(ai + 1, args.alpha_parts))
                for hj in range(args.h_parts):
                    hh = I(Q(hj, args.K_min * args.h_parts), Q(hj + 1, args.K_min * args.h_parts))
                    for eta in (I(Q(-1, 2), Q(0)), I(Q(0), Q(1, 2))):
                        visit(branch, r, aa, hh, eta)
            print(branch, r, 'boxes', count, 'failures', len(failures), flush=True)
    report = {'status': 'ALL_BOXES_POSITIVE' if not failures else 'OPEN',
              'alpha': list(map(str, ALPHA)), 'K_min': args.K_min,
              'branch_alpha': {k: list(map(str, v)) for k, v in BRANCH_ALPHA.items()},
              'strict_lower_bounds': {k: str(v) for k, v in FLOORS.items()},
              'alpha_parts': args.alpha_parts, 'h_parts': args.h_parts,
              'box_count': count, 'minima': minima, 'failure_count': len(failures),
              'terminal_box_count': len(leaves), 'terminal_boxes_sha256': hashlib.sha256(json.dumps(leaves).encode()).hexdigest(),
              'first_failures': failures[:20], 'elapsed_seconds': time.monotonic() - started,
              'source_sha256': {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [Path(__file__), ROOT/'tools/audit_v04_uniform_tail.py']}}
    args.report.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(report['status'], count, {k: v['approximation'] for k, v in minima.items()}, flush=True)
    return 0 if not failures else 1


if __name__ == '__main__':
    raise SystemExit(main())
