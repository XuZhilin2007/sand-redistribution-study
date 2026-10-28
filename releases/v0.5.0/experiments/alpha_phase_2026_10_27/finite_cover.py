"""Exact all-time finite complement on [123/500,7/24], excluding K=6,10.

Closed cells cover every s1 envelope segment, including both branch values
at selector switches. Future boundary words are over-approximated by a binary
tree. Impossible paths alone are pruned; each remaining loss is certified.
An independent square-completion tail covers all subsequent integer times.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
import audit_v04_interval_v3 as P
from phase_probe import model
from verify_phase import polynomial_states, discrepancies, tail_polynomials

DOMAIN = (Q(123, 500), Q(7, 24))
RHO = (Q(1), Q(-1))


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def positive(p, lo, hi, strict=True):
    return P.certify_sign(p, lo, hi, strict=strict, max_depth=16)[0]


def envelope(k, domain=DOMAIN):
    """Upper hull of exact affine discrepancies, with a tie coverage check."""
    states, weights, h = polynomial_states(k, [k // 2])
    ds = discrepancies(states[1], k)
    sorted_lines = sorted((p[1] if len(p) > 1 else Q(0), p[0] if p else Q(0), j)
                          for j, p in enumerate(ds, 1))
    unique = []
    for slope, intercept, index in sorted_lines:
        line = (slope, intercept, index)
        if unique and unique[-1][0] == slope:
            if (intercept, -index) > (unique[-1][1], -unique[-1][2]):
                unique[-1] = line
        else:
            unique.append(line)
    hull = []
    for line in unique:
        start = None
        while hull:
            last, old_start = hull[-1]
            start = (last[1] - line[1]) / (line[0] - last[0])
            if old_start is None or start > old_start:
                break
            hull.pop()
        if not hull:
            start = None
        hull.append((line, start))
    cells = []
    for i, (line, start) in enumerate(hull):
        lo = max(domain[0], start) if start is not None else domain[0]
        finish = hull[i + 1][1] if i + 1 < len(hull) else None
        hi = min(domain[1], finish) if finish is not None else domain[1]
        if lo < hi:
            cells.append((lo, hi, line[2]))
    require(cells[0][0] == domain[0] and cells[-1][1] == domain[1], 'envelope endpoints')
    for i, (lo, hi, s) in enumerate(cells):
        require(i == 0 or lo == cells[i - 1][1], 'envelope has gap')
        for a in (lo, hi):
            actual = 1 + max(range(k), key=lambda j: (P.value(ds[j], a), -j))
            allowed = {s}
            if i and a == lo:
                allowed.add(cells[i - 1][2])
            if i + 1 < len(cells) and a == hi:
                allowed.add(cells[i + 1][2])
            require(actual in allowed, 'uncovered isolated minimum-index tie')
    return cells


@lru_cache(maxsize=None)
def branch_data(k, s1):
    states, weights, h = polynomial_states(k, [k // 2, s1])
    u0, u1, u2 = [P.full_loss(p, k) for p in states]
    ds1, ds2 = discrepancies(states[1], k), discrepancies(states[2], k)
    h2, h1, h0 = tail_polynomials(k, states, weights)
    return states, weights, h, (u0, u1, u2), ds1, ds2, (h2, h1, h0)


def finite_tree(k, states, weights, lo, hi, last_time):
    """All possibly legal boundary words through last_time, with exact ties."""
    h, m = Q(1, k), (5 * k + 5) // 6
    u2 = P.full_loss(states[2], k)
    obligations = []
    def visit(state, t, path):
        if t >= 3:
            gap = P.sub(P.full_loss(state, k), u2)
            if not positive(gap, lo, hi):
                return False
            obligations.append((f'U{t}>U2:{path}', gap, True))
        if t == last_time:
            return True
        x_minus_h = P.sub(state[m - 1], (h,))
        for s in (m - 1, m):
            # x=h belongs to ell. A weak upper bound x<=h discards m;
            # discarding ell requires a strictly positive lower bound.
            impossible = P.scale(x_minus_h, -1) if s == m else x_minus_h
            strict = s == m - 1
            if positive(impossible, lo, hi, strict):
                obligations.append((f'prune:{path},{s}', impossible, strict))
            elif not visit(P.update_full(state, s, weights), t + 1, path + (s,)):
                return False
        return True
    return visit(states[2], 2, ()), obligations


def try_cell(k, s1, lo, hi, first_tail_fixed=None):
    states, weights, h, uu, ds1, ds2, hs = branch_data(k, s1)
    m = (5 * k + 5) // 6
    obligations = []
    def check(p, label, strict=True):
        if not positive(p, lo, hi, strict):
            return False
        obligations.append((label, p, strict))
        return True
    # Nonstrict global license permits both limiting branches at each end.
    # The envelope check proves the actual tie branch is among covered cells.
    for j, d in enumerate(ds1, 1):
        if j != s1 and not check(P.sub(ds1[s1 - 1], d), f's1 versus {j}', False):
            return None
    if s1 not in (m - 1, m):
        # Every internal coordinate must be strictly dominated by at least
        # one boundary coordinate. Theorem A excludes all j>m separately.
        for j in range(m - 2):
            candidates = [P.sub(ds2[s - 1], ds2[j]) for s in (m - 1, m)]
            chosen = next((p for p in candidates if positive(p, lo, hi)), None)
            if chosen is None:
                return None
            obligations.append((f'entry versus {j+1}', chosen, True))
    for j in (0, 1):
        if not check(P.sub(uu[j], uu[2]), f'U{j}>U2'):
            return None
    h2, h1, h0 = hs
    for first_tail in ([first_tail_fixed] if first_tail_fixed is not None else range(3, 9)):
        w = P.power(RHO, first_tail - 2)
        hp = P.scale(P.add(P.scale(P.mul(h2, w), 2), h1), -1)
        hw = P.total((P.mul(h2, P.mul(w, w)), P.mul(h1, w), h0))
        if not (positive(hp, lo, hi) and positive(hw, lo, hi)):
            continue
        if first_tail > 3:
            good, tree_obligations = finite_tree(k, states, weights, lo, hi, first_tail - 1)
            if not good:
                continue
        else:
            tree_obligations = []
        obligations.extend(tree_obligations + [('negative tail derivative', hp, True), ('positive tail', hw, True)])
        digest = hashlib.sha256(json.dumps([(label, strict, list(map(str, p)))
                                          for label, p, strict in obligations]).encode()).hexdigest()
        return {'K': k, 's1': s1, 'lo': str(lo), 'hi': str(hi), 'mode': 'closed',
                'first_tail_time': first_tail, 'obligations': len(obligations),
                'obligations_sha256': digest}
    return None


def cover(k, s1, lo, hi, depth=0):
    result = try_cell(k, s1, lo, hi)
    if result is not None:
        return [result]
    require(depth < 12, f'unresolved K={k},s1={s1},interval={lo},{hi}')
    mid = (lo + hi) / 2
    return cover(k, s1, lo, mid, depth + 1) + cover(k, s1, mid, hi, depth + 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-k', type=int, default=203)
    parser.add_argument('--report', type=Path, default=Path(__file__).with_name('finite_cover.json'))
    args = parser.parse_args()
    started = time.monotonic()
    records, targets = [], [k for k in range(7, args.max_k + 1) if k != 10]
    for k in targets:
        rows = []
        for lo, hi, s1 in envelope(k):
            rows.extend(cover(k, s1, lo, hi))
        require(Q(rows[0]['lo']) == DOMAIN[0] and Q(rows[-1]['hi']) == DOMAIN[1], 'coverage endpoint')
        require(all(Q(a['hi']) == Q(b['lo']) for a, b in zip(rows, rows[1:])), 'coverage gap')
        records.extend(rows)
        branch_data.cache_clear()
        if k <= 20 or k % 10 == 0:
            print('K', k, 'cells', len(rows), 'tail<=', max(r['first_tail_time'] for r in rows), flush=True)
    extra10 = []
    for lo, hi, s1 in envelope(10, (Q(1, 4), DOMAIN[1])):
        extra10.extend(cover(10, s1, lo, hi))
    report = {'status': 'VERIFIED_FINITE_ALL_TIME_COVER', 'alpha': list(map(str, DOMAIN)),
              'targets': targets, 'records': records, 'cell_count': len(records),
              'K10_upper_domain': ['1/4', str(DOMAIN[1])], 'K10_upper_records': extra10,
              'obligation_count': sum(r['obligations'] for r in records),
              'elapsed_seconds': time.monotonic() - started,
              'source_sha256': {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [Path(__file__), Path(__file__).with_name('verify_phase.py'),
                                          ROOT / 'tools/audit_v04_interval_v3.py']}}
    args.report.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('VERIFIED', len(targets), 'K', len(records), 'cells', report['elapsed_seconds'], 'seconds', flush=True)


if __name__ == '__main__':
    main()
