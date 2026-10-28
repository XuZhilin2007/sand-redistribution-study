"""Exact post-release phase results. Never writes frozen sources or evidence.

Uses two existing, different sign methods (Sturm and rational Bernstein),
reconstructs all trajectory polynomials, and licenses the minimum-index rule.
The analytic all-time implications are in docs/alpha_phase_2026-10-27/RESULTS.md.
No sample, status string, or finite replay licenses an infinite quantifier.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from functools import reduce
import hashlib
import json
from math import gcd, lcm
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'tools'), str(ROOT / 'experiments/v04_w_interval')]
import audit_v04_interval_v3 as P
import symkernel as S
from phase_probe import model, replay, select

RHO = (Q(1), Q(-1))
OBLIGATIONS = []


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def sign(poly, lo, hi, label, strict=True, left_open=False):
    """Every accepted sign is checked by BOTH exact algorithms."""
    poly = P.trim(poly)
    sturm = ((S.pos_leftopen if left_open else S.pos_closed) if strict else
             (S.nonneg_leftopen if left_open else S.nonneg_closed))
    # The audited kernel deliberately treats the identically zero polynomial
    # as OPEN; zero is nevertheless a direct nonnegative identity.
    a = (not poly and not strict) or sturm(poly, lo, hi)
    b, leaves, depth = P.certify_sign(poly, lo, hi, strict, left_open, max_depth=20)
    require(a and b, label)
    OBLIGATIONS.append({'label': label, 'lo': str(lo), 'hi': str(hi),
                        'strict': strict, 'left_open': left_open,
                        'coefficients_ascending': list(map(str, poly)),
                        'bernstein_leaves': leaves, 'bernstein_depth': depth})


def can_sign(poly, lo, hi, strict=True, left_open=False):
    poly = P.trim(poly)
    f = ((S.pos_leftopen if left_open else S.pos_closed) if strict else
         (S.nonneg_leftopen if left_open else S.nonneg_closed))
    return (not poly and not strict) or f(poly, lo, hi)


def polynomial_states(k, word):
    p, weights, h = model(k, True)
    states = [[(x,) for x in p]]
    for s in word:
        states.append(P.update_full(states[-1], s, weights))
    return states, weights, h


def discrepancies(state, k):
    result, acc = [], ()
    for j, poly in enumerate(state, 1):
        acc = P.add(acc, poly)
        result.append(P.sub(acc, (Q(j, k),)))
    return result


def license_selector(state, s, lo, hi, label, left_open=False):
    ds = discrepancies(state, len(state))
    for j, other in enumerate(ds, 1):
        if j != s:
            sign(P.sub(ds[s - 1], other), lo, hi, label + f':versus-{j}',
                 strict=j < s, left_open=left_open)


def losses(states):
    return [P.full_loss(p, len(p)) for p in states]


def primitive(poly):
    poly = P.trim(poly)
    common = lcm(*(x.denominator for x in poly))
    ints = [int(x * common) for x in poly]
    divisor = reduce(gcd, ints)
    ints = [x // divisor for x in ints]
    if ints[-1] < 0:
        ints = [-x for x in ints]
    return ints


def irreducibility_mod_prime(coefficients):
    """For degrees 3 or 5, exclude all possible small-degree factors mod p."""
    degree = len(coefficients) - 1
    if degree not in (3, 5):
        return None
    for prime in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        if coefficients[-1] % prime == 0:
            continue
        if any(sum(c * pow(x, j, prime) for j, c in enumerate(coefficients)) % prime == 0
               for x in range(prime)):
            continue
        if degree == 5:
            quadratic_factor = False
            for a in range(prime):
                for b in range(prime):
                    rem = [v % prime for v in coefficients]
                    for j in range(degree, 1, -1):
                        factor = rem[j]
                        rem[j - 1] = (rem[j - 1] - factor * a) % prime
                        rem[j - 2] = (rem[j - 2] - factor * b) % prime
                        rem[j] = 0
                    if rem[0] == rem[1] == 0:
                        quadratic_factor = True
                        break
                if quadratic_factor:
                    break
            if quadratic_factor:
                continue
        return prime
    return None


def root_record(poly, lo, hi, label):
    """Unique simple root: Sturm count plus a Bernstein derivative sign."""
    poly = P.trim(poly)
    require(S.count_roots_open(poly, lo, hi) == 1, label + ':unique root')
    vl, vh = P.value(poly, lo), P.value(poly, hi)
    require(vl * vh < 0, label + ':opposite signs')
    derivative = S.pderiv(poly)
    direction = 1 if vh > vl else -1
    sign(P.scale(derivative, direction), lo, hi, label + ':strict monotonicity')
    aa, bb = lo, hi
    for _ in range(100):
        mid = (aa + bb) / 2
        vm = P.value(poly, mid)
        if vm == 0:
            aa = bb = mid
            break
        if (vm > 0) == (vl > 0):
            aa = mid
        else:
            bb = mid
    # Remove powers of alpha only for reporting the defining polynomial.
    reduced = list(poly)
    while reduced and not reduced[0]:
        reduced.pop(0)
    coefficients = primitive(reduced)
    return {'label': label, 'unique_in': [str(lo), str(hi)],
            'bracket': [str(aa), str(bb)], 'approximation': float((aa + bb) / 2),
            'primitive_polynomial_ascending': coefficients,
            'irreducible_mod_prime': irreducibility_mod_prime(coefficients),
            'loss_polynomial_ascending': list(map(str, poly))}


def exclude_by_U3(k, s1, lo, hi, depth=0):
    """Prove U3<U2 on (lo,hi] with all feasible s2, without guessing a word.

    A candidate s2 is discarded only when an earlier index weakly dominates
    it throughout the cell, or another index strictly dominates it. Other
    candidates must have strictly negative U3-U2. Rational splitting covers
    selector switches and does not assume cross-branch continuity.
    """
    states, weights, h = polynomial_states(k, [k // 2, s1])
    license_selector(states[1], s1, lo, hi, f'K{k}:s1={s1}', True)
    ds = discrepancies(states[2], k)
    u2 = P.full_loss(states[2], k)
    pending = []
    obligations = []
    for s2 in range(1, k + 1):
        impossible = next((j for j in range(1, k + 1) if j != s2 and
                           can_sign(P.sub(ds[j - 1], ds[s2 - 1]), lo, hi,
                                    strict=j > s2, left_open=True)), None)
        if impossible is not None:
            obligations.append((P.sub(ds[impossible - 1], ds[s2 - 1]),
                                impossible > s2, f'exclude-s2={s2}-by-{impossible}'))
        else:
            p3 = P.update_full(states[2], s2, weights)
            improvement = P.sub(u2, P.full_loss(p3, k))
            if can_sign(improvement, lo, hi, left_open=True):
                obligations.append((improvement, True, f'U2>U3:s2={s2}'))
            else:
                pending.append(s2)
    if pending:
        require(depth < 14, f'K{k} exclusion unresolved {lo},{hi}: {pending}')
        mid = (lo + hi) / 2
        return (exclude_by_U3(k, s1, lo, mid, depth + 1) +
                exclude_by_U3(k, s1, mid, hi, depth + 1))
    for poly, strict, label in obligations:
        sign(poly, lo, hi, f'K{k}:{label}', strict=strict, left_open=True)
    return [{'K': k, 's1': s1, 'domain': [str(lo), str(hi)], 'mode': 'left-open',
             'counterexample_time': 3, 'feasible_s2': [int(x[2].split('=')[-1])
                                                      for x in obligations if x[2].startswith('U2>')]}]


def tail_polynomials(k, states, weights):
    """Full-vector square completion anchored at p2; E1=0 checked exactly."""
    h = Q(1, k)
    m = (5 * k + 5) // 6
    b = m * h
    c, q = b * (2 - b), weights[m - 1]
    gamma = P.scale(P.sub((c,), states[2][m - 1]), 1 / (1 - q))
    error = [P.sub(states[2][i], P.scale(gamma, weights[i])) for i in range(m - 1)]
    require(not P.total(error), 'zero residual mass')
    base = [c * y / (1 - q) for y in weights[:m - 1]] + [Q(0)]
    base += [x[0] for x in states[0][m:]]
    direction = [-y / (1 - q) for y in weights[:m - 1]] + [Q(1)] + [Q(0)] * (k - m)
    error += [()] * (k - m + 1)
    a = k * sum(x * x for x in direction)
    bb = 2 * k * sum((x - h) * y for x, y in zip(base, direction))
    cc = k * sum((x - h) ** 2 for x in base)
    dd = P.scale(P.total(P.scale(e, v) for e, v in zip(error, direction)), k)
    ff = P.scale(P.total(P.scale(e, v - h) for e, v in zip(error, base)), k)
    ee = P.scale(P.total(P.mul(e, e) for e in error), k)
    h2 = P.sub(ee, P.scale(P.mul(dd, dd), 1 / a))
    h1 = P.sub(P.scale(ff, 2), P.scale(dd, bb / a))
    h0 = P.sub((cc - bb * bb / (4 * a),), P.full_loss(states[2], k))
    return h2, h1, h0


def certify_tail(k, states, weights, lo, hi, first_tail_time, label):
    h2, h1, h0 = tail_polynomials(k, states, weights)
    w = P.power(RHO, first_tail_time - 2)
    # Gram gives H2>=0 on actual states. These two signs cover all z<=w.
    hp = P.add(P.scale(P.mul(h2, w), 2), h1)
    hw = P.total((P.mul(h2, P.mul(w, w)), P.mul(h1, w), h0))
    sign(P.scale(hp, -1), lo, hi, label + ':Hprime<0')
    sign(hw, lo, hi, label + ':H(w)>0')


def central_results():
    cells = []
    for k, pieces in ((6, [(3, Q(0), Q(4, 27)), (4, Q(4, 27), Q(19, 100))]),
                      (10, [(5, Q(0), Q(4, 45)), (6, Q(4, 45), Q(6, 25))])):
        for s1, lo, hi in pieces:
            cells.extend(exclude_by_U3(k, s1, lo, hi))

    # K6 central branch. Once x2>h, cq>h makes selector m permanent.
    states6, g6, h6 = polynomial_states(6, [3, 4, 5])
    u6 = losses(states6)
    lo, hi = Q(19, 100), Q(1, 3)
    license_selector(states6[1], 4, lo, hi, 'K6 central s1')
    license_selector(states6[2], 5, lo, hi, 'K6 central s2')
    require(Q(35, 36) * g6[4] > h6, 'K6 permanent m fixed point above h')
    sign(P.sub(u6[0], u6[2]), lo, hi, 'K6 central U0>U2')
    lower6 = root_record(P.sub(u6[3], u6[2]), Q(19, 100), Q(1, 5), 'a6')
    sign(P.sub(u6[3], u6[2]), Q(1, 5), hi, 'K6 U3>U2 past a6 bracket')
    upper = root_record(P.sub(u6[1], u6[2]), lo, hi, 'b6')

    # High alpha K6 exclusion: the exact s1 changes at 1/3.
    high6, gg, hh = polynomial_states(6, [3, 5])
    uh = losses(high6)
    license_selector(high6[1], 5, Q(1, 3), Q(1), 'K6 high s1', True)
    sign(P.sub(uh[2], uh[1]), Q(1, 3), Q(2, 3), 'K6 high U1<U2')
    sign(P.sub(uh[2], uh[0]), Q(2, 3), Q(1), 'K6 high U0<U2')

    # K10 true word on (6/25,1/4], with a loss crossing, not a selector switch.
    s10, g10, h10 = polynomial_states(10, [5, 7, 8])
    u10 = losses(s10)
    license_selector(s10[1], 7, Q(6, 25), Q(1, 4), 'K10 central s1', True)
    license_selector(s10[2], 8, Q(6, 25), Q(1, 4), 'K10 central s2')
    lower = root_record(P.sub(u10[3], u10[2]), Q(6, 25), Q(1, 4), 'a10')
    for j in (0, 1):
        sign(P.sub(u10[j], u10[2]), Q(6, 25), Q(1, 4), f'K10 central U{j}>U2')
    certify_tail(10, s10, g10, Q(6, 25), Q(1, 4), 4, 'K10 central t>=4')
    return {'roots': {'a6': lower6, 'a10': lower, 'b6': upper},
            'low_alpha_exclusion_cells': cells,
            'theorems': ['W6=[a6,b6]', 'W10 intersect (6/25,1/4]=[a10,1/4]',
                         'intersection of all W_K is contained in [a10,b6]']}


def disconnected_results():
    # A distant K10 component, certified directly for all t>=3.
    states, weights, h = polynomial_states(10, [5, 8])
    lo, hi = Q(21, 25), Q(43, 50)
    license_selector(states[1], 8, lo, hi, 'K10 distant s1=ell')
    # Selection s1=ell is already boundary entry; Theorem A gives residence.
    uu = losses(states)
    for j in (0, 1):
        sign(P.sub(uu[j], uu[2]), lo, hi, f'K10 distant U{j}>U2')
    certify_tail(10, states, weights, lo, hi, 3, 'K10 distant t>=3')

    # K13 immediately to the right of an exact s1 selector discontinuity.
    s13, g13, h13 = polynomial_states(13, [6, 10, 11])
    lo13, hi13 = Q(13, 40), Q(1, 3)
    license_selector(s13[1], 10, lo13, hi13, 'K13 right branch s1', True)
    license_selector(s13[2], 11, lo13, hi13, 'K13 right branch s2')
    m = 11
    c = Q(m, 13) * (2 - Q(m, 13))
    require(c * g13[m - 1] > h13, 'K13 permanent m fixed point above h')
    u13 = losses(s13)
    for j in (0, 1, 3):
        sign(P.sub(u13[j], u13[2]), lo13, hi13, f'K13 right branch U{j}>U2')

    actual, word, uu13 = replay(13, lo13, 3)
    require(word[:2] == [6, 9], 'K13 tie is resolved to 9')
    require(uu13[1] < uu13[2], 'K13 tie fails')
    right_limit_u2 = P.value(u13[2], lo13)
    jump = right_limit_u2 - uu13[2]
    return {'K10_certified_distant_interval': ['21/25', '43/50'],
            'K13_certified_half_open_interval': ['13/40', '1/3'],
            'K13_tie_word': word, 'K13_tie_U1_minus_U2': str(uu13[1] - uu13[2]),
            'K13_right_limit_U2_minus_actual_U2': str(jump),
            'theorems': ['W10 is disconnected', 'W13 is disconnected and is not closed in (0,1)']}


def exact_witnesses():
    rows = []
    for k, alpha, time_index in ((9, Q(1, 5), 3), (6, Q(7, 24), 1),
                                  (10, Q(1, 2), 1), (10, Q(9, 10), 4),
                                  (13, Q(81, 250), 1), (13, Q(13, 40), 1)):
        states, word, uu = replay(k, alpha, max(4, time_index))
        require(uu[time_index] < uu[2], f'exact witness {k},{alpha}')
        if k == 10 and alpha == Q(9, 10):
            require(all(uu[j] > uu[2] for j in (0, 1, 3)), 'fourth-time witness passes earlier comparisons')
        rows.append({'K': k, 'alpha': str(alpha), 'word': word,
                     'U': list(map(str, uu)), 'witness_time': time_index,
                     'gap': str(uu[time_index] - uu[2])})

    # Fraction-vector replay is independent of polynomial expansion. Check
    # the global selectors and all vector entries for several branch words.
    checks = 0
    for k in (6, 7, 8, 9, 10, 11, 13, 16, 32):
        for alpha in (Q(19, 100), Q(6, 25), Q(1, 4), Q(13, 40), Q(17, 20)):
            exact, word, uu = replay(k, alpha, 5)
            polys, gg, hh = polynomial_states(k, word)
            for v1, v2 in zip(exact, polys):
                require(v1 == [P.value(poly, alpha) for poly in v2], 'independent vector polynomial replay')
            checks += 1
    require(rows[0]['gap'] == '-83303117/5978711250', 'published left witness')
    require(rows[1]['gap'] == '-161/11943936', 'published right witness')
    # These controls reject two tempting incorrect statements.
    require(not can_sign((Q(1), Q(-5)), Q(1, 4), Q(1, 3)), 'negative sign control')
    return {'witnesses': rows, 'independent_polynomial_vector_cases': checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, default=Path(__file__).with_name('exact_results.json'))
    args = parser.parse_args()
    started = time.monotonic()
    central = central_results()
    print('central results and global exclusion: PASS', flush=True)
    disconnected = disconnected_results()
    print('disconnected and nonclosed fixed-K sets: PASS', flush=True)
    witnesses = exact_witnesses()
    sources = [Path(__file__), Path(__file__).with_name('phase_probe.py'),
               ROOT / 'tools/audit_v04_interval_v3.py',
               ROOT / 'experiments/v04_w_interval/symkernel.py',
               ROOT / 'docs/XO_ALPHA_EXTENSION_PROOF.md']
    report = {'status': 'VERIFIED_POST_RELEASE_PHASE_RESULTS',
              'scope': 'Only the exact theorems explicitly listed; universal inner interval inherited from v0.4.0',
              'central': central, 'disconnected': disconnected, **witnesses,
              'sign_obligation_count': len(OBLIGATIONS), 'sign_obligations': OBLIGATIONS,
              'elapsed_seconds': time.monotonic() - started,
              'source_sha256': {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sources}}
    args.report.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('exact signs', len(OBLIGATIONS), 'seconds', report['elapsed_seconds'])
    for name, r in central['roots'].items():
        print(name, r['approximation'], 'primitive ascending', r['primitive_polynomial_ascending'])
    print('wrote', args.report)


if __name__ == '__main__':
    main()
