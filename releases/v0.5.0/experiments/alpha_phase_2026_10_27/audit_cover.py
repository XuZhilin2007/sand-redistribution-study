"""Read the finite covering witness and recheck it with exact Sturm signs.

The producer used Bernstein coefficients. This pass checks the recorded tail
time, reconstructs the proof from the model, checks coverage and the exact
tie winners, and exercises rejection controls. Shared trajectory algebra and
analytic residence remain explicit common premises.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import time

import finite_cover as F
from verify_phase import S, P, polynomial_states, discrepancies


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def sturm_positive(poly, lo, hi, strict=True):
    poly = P.trim(poly)
    if not poly:
        return not strict
    return (S.pos_closed if strict else S.nonneg_closed)(poly, lo, hi)


def check_coverage(rows, k, lo, hi):
    require(rows and all(r['K'] == k for r in rows), 'empty/wrong K covering')
    ordered = sorted(rows, key=lambda r: (Q(r['lo']), Q(r['hi'])))
    require(Q(ordered[0]['lo']) == lo and Q(ordered[-1]['hi']) == hi, 'wrong endpoints')
    require(all(Q(a['hi']) == Q(b['lo']) for a, b in zip(ordered, ordered[1:])), 'gap or overlap')
    require(all(r['mode'] == 'closed' and Q(r['lo']) < Q(r['hi']) for r in ordered), 'invalid cell')
    # Actual minimum-index winner at every cell boundary must be among the
    # cells that include that point. This also handles a three-way tie.
    state1 = polynomial_states(k, [k//2])[0][1]
    ds = discrepancies(state1, k)
    for alpha in {Q(r[key]) for r in ordered for key in ('lo', 'hi')}:
        actual = 1 + max(range(k), key=lambda j: (P.value(ds[j], alpha), -j))
        covering = [r for r in ordered if Q(r['lo']) <= alpha <= Q(r['hi'])]
        require(any(r['s1'] == actual for r in covering), f'uncertified exact tie K={k},alpha={alpha}')
    return ordered


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path(__file__).with_name('finite_cover.json'))
    parser.add_argument('--report', type=Path, default=Path(__file__).with_name('finite_audit.json'))
    args = parser.parse_args()
    started = time.monotonic()
    witness = json.loads(args.input.read_text(encoding='utf-8'))
    targets = set(range(7, 204)) - {10}
    require(witness['alpha'] == ['123/500', '7/24'], 'wrong finite interval')
    require(set(witness['targets']) == targets and len(witness['targets']) == len(targets), 'wrong target set')
    grouped = defaultdict(list)
    for row in witness['records']:
        grouped[row['K']].append(row)
    require(set(grouped) == targets, 'missing or extra K')
    F.positive = sturm_positive
    count = 0
    for k in sorted(targets | {10}):
        rows = grouped[k] if k != 10 else witness['K10_upper_records']
        lo = Q(123, 500) if k != 10 else Q(1, 4)
        rows = check_coverage(rows, k, lo, Q(7, 24))
        for row in rows:
            require(type(row['first_tail_time']) is int and 3 <= row['first_tail_time'] <= 8, 'invalid tail time')
            rebuilt = F.try_cell(k, row['s1'], Q(row['lo']), Q(row['hi']), row['first_tail_time'])
            require(rebuilt is not None, f'not certified K={k}: {row}')
            require(rebuilt['first_tail_time'] == row['first_tail_time'], 'wrong time coverage')
            # Exact strictness inside the cell, even if endpoints tie.
            ds = discrepancies(polynomial_states(k, [k//2])[0][1], k)
            for j in range(row['s1']-1):
                require(P.trim(P.sub(ds[row['s1']-1], ds[j])), 'permanently tied smaller index')
            count += 1
        F.branch_data.cache_clear()
        if k % 40 == 0:
            print('Sturm audit through K', k, 'cells', count, flush=True)

    # End-to-end controls corresponding to possible proof mistakes.
    require(F.try_cell(6, 4, Q(299, 1000), Q(3, 10)) is None, 'reject omitted time-1 obligation')
    require(F.try_cell(10, 7, Q(241, 1000), Q(245, 1000)) is None, 'reject true U3 counterexample')
    require(F.try_cell(9, 6, Q(123, 500), Q(1, 4), 3) is None, 'reject tail starting too early')
    bad32 = [dict(row) for row in grouped[32]]
    for row in bad32:
        if Q(row['lo']) <= Q(1, 4) <= Q(row['hi']):
            row['s1'] = 22
    try:
        check_coverage(bad32, 32, Q(123, 500), Q(7, 24))
    except AssertionError:
        pass
    else:
        raise AssertionError('missing smallest-index triple tie was accepted')
    report = {'status': 'VERIFIED_FINITE_COVER_BY_STURM', 'verified_targets': sorted(targets),
              'K10_upper_verified': True, 'verified_cells': count,
              'rejection_controls': 4, 'elapsed_seconds': time.monotonic()-started,
              'witness_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
              'source_sha256': {str(p.relative_to(F.ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [Path(__file__), Path(F.__file__), Path(__file__).with_name('verify_phase.py'),
                                          F.ROOT/'experiments/v04_w_interval/symkernel.py']}}
    args.report.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('PASS', count, 'cells; 4 rejection controls;', report['elapsed_seconds'], 'seconds', flush=True)


if __name__ == '__main__':
    main()
