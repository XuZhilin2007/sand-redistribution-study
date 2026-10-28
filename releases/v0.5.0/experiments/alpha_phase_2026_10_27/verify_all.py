"""Reproduce and assemble the post-release sharp universal alpha theorem.

Runs the actual proof checks, checks explicit target coverage and root order,
and verifies that v0.4.0 and the prior review/corrections remain unchanged.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((HERE/name).read_text(encoding='utf-8'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-baseline-replay', action='store_true')
    parser.add_argument('--report', type=Path, default=HERE/'verification.json')
    args = parser.parse_args()
    started = time.monotonic()
    sys.path.insert(0, str(ROOT/'tools'))
    import verify_v040 as baseline
    expected = baseline.manifest_files()
    require(not baseline.check_hashes(expected), 'frozen baseline initial hashes')
    preserved_paths = [ROOT/'docs/V040_MANIFEST.json', ROOT/'docs/READING_NOTES.md', ROOT/'docs/MODEL_ASSUMPTIONS.md']
    preserved_paths += [p for sub in ('review_2026-10-27', 'corrections_2026-10-27')
                        for p in (ROOT/'docs'/sub).rglob('*') if p.is_file()]
    preserved = {str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in preserved_paths}
    checks = []
    commands = [
        [str(HERE/'verify_phase.py')],
        [str(HERE/'audit_phase.py')],
        [str(HERE/'finite_cover.py'), '--max-k', '203'],
        [str(HERE/'audit_cover.py')],
        [str(HERE/'uniform_probe.py'), '--alpha-parts', '32', '--h-parts', '4', '--K-min', '204'],
    ]
    if not args.skip_baseline_replay:
        commands.append([str(ROOT/'tools/verify_v040.py'), '--report', str(HERE/'baseline_verification.json')])
    for command in commands:
        print('Checking', Path(command[0]).name, flush=True)
        tick = time.monotonic()
        completed = subprocess.run([sys.executable, '-B', *command], cwd=ROOT,
                                   capture_output=True, encoding='utf-8')
        check = {'command': [sys.executable, '-B', *command], 'exit_code': completed.returncode,
                 'seconds': time.monotonic()-tick, 'stdout': completed.stdout, 'stderr': completed.stderr}
        checks.append(check)
        if completed.returncode:
            print(completed.stdout, completed.stderr, flush=True)
            raise AssertionError('failed check ' + command[0])
        print('  PASS', round(check['seconds'], 3), 'seconds', flush=True)

    exact, algebra, finite, second, boxes = [read(name) for name in
        ('exact_results.json', 'algebra_audit.json', 'finite_cover.json', 'finite_audit.json', 'uniform_boxes.json')]
    for record in (exact, algebra, finite, second, boxes):
        for name, expected_digest in record['source_sha256'].items():
            require(digest(ROOT/name) == expected_digest, 'source binding ' + name)
    roots = exact['central']['roots']
    a6, aa, bb = [[Q(v) for v in roots[name]['bracket']] for name in ('a6', 'a10', 'b6')]
    require(a6[1] < Q(123, 500) < aa[0] < aa[1] < Q(1, 4) < bb[0] < bb[1] < Q(7, 24), 'exact root order')
    require(roots['a10']['irreducible_mod_prime'] == 29 and roots['b6']['irreducible_mod_prime'] == 11,
            'boundary degrees 5 and 3')
    targets = set(range(7, 204))-{10}
    require(set(finite['targets']) == set(second['verified_targets']) == targets, 'complete finite complement')
    require(second['K10_upper_verified'] and second['witness_sha256'] == digest(HERE/'finite_cover.json'), 'independent finite audit binding')
    require(algebra['formal_identity_count'] == 60 and len(algebra['exact_vector_cases']) == 95, 'algebra obligations')
    require(boxes['status'] == 'ALL_BOXES_POSITIVE' and boxes['failure_count'] == 0, 'uniform sign proof')
    require(boxes['K_min'] == 204 and boxes['alpha'] == ['123/500', '7/24'], 'uniform domain and no K gap')
    require(boxes['branch_alpha'] == {'M': ['123/500', '13/50'], 'U': ['1/4', '7/24']}, 'true branch coverage')
    for name, floor in boxes['strict_lower_bounds'].items():
        require(Q(boxes['minima'][name]['lower']) > Q(floor), 'uniform positive margin ' + name)
    require(not baseline.check_hashes(expected), 'frozen baseline final hashes')
    require(all(digest(ROOT/name) == expected_digest for name, expected_digest in preserved.items()), 'review/corrections unchanged')
    reports = ('exact_results.json', 'algebra_audit.json', 'finite_cover.json', 'finite_audit.json', 'uniform_boxes.json')
    report = {'status': 'VERIFIED_SHARP_UNIVERSAL_ALPHA_INTERVAL',
              'claim': 'For the frozen canonical XO model, intersection_{K>=6} W_K = W_6 intersect W_10 = [a10,b6].',
              'boundary_roots': {name: roots[name] for name in ('a10', 'b6')},
              'fixed_K_results': exact['central']['theorems'][:2] + exact['disconnected']['theorems'],
              'all_other_K': 'For every K>=7 with K!=10: [123/500,7/24] is contained in W_K.',
              'finite_range': [7, 203], 'uniform_range': [204, 'infinity'],
              'baseline_replayed': not args.skip_baseline_replay, 'baseline_files_unchanged': len(expected),
              'prior_review_and_corrections_unchanged': preserved,
              'checks': checks, 'reports_sha256': {name: digest(HERE/name) for name in reports},
              'source_sha256': {str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in sorted(HERE.glob('*.py'))},
              'elapsed_seconds': time.monotonic()-started,
              'qualification': 'Computer-assisted mathematical proof under the documented analytic premises; no claim of human peer review or proof-assistant verification.'}
    args.report.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(report['status'], 'seconds', report['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
