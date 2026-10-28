"""Focused evidence-integrity rejection tests; never run proof generation.

All intentionally changed inputs live in NEW test copies, outside the frozen
package. The tested package and original historical evidence stay unchanged.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

from v050_common import (REPORTS, canonical_digest, digest, mathematical_payload,
    new_output_directory, read_json, require, safe_extract, verify_inventory, write_new_json)
from reproduce_v050 import preflight, validate_claims


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    package = args.package.resolve()
    manifest, pin = verify_inventory(package)
    output = new_output_directory(package, args.out_dir)
    tests = []
    def check(name, condition):
        require(condition, 'test failed: '+name)
        tests.append({'test': name, 'status': 'PASS'})
    def refused(name, operation):
        try:
            operation()
        except (ValueError, OSError, KeyError):
            tests.append({'test': name, 'status': 'PASS'})
        else:
            raise AssertionError('dangerous action accepted: '+name)

    preflight(package, pin)
    check('complete unmodified package accepted', True)
    refused('wrong externally pinned manifest rejected', lambda: preflight(package, '0'*64))
    refused('output inside package rejected', lambda: new_output_directory(package, package/'generated'))
    refused('normalized output path inside package rejected', lambda: new_output_directory(package, package/'docs'/'..'/'generated'))
    existing = output/'existing-evidence'
    existing.mkdir()
    (existing/'sentinel.txt').write_bytes(b'keep evidence bytes\n')
    before = digest(existing/'sentinel.txt')
    refused('existing output rejected', lambda: new_output_directory(package, existing))
    check('existing evidence unchanged', digest(existing/'sentinel.txt') == before)

    for label, relative, action in (
        ('changed proof rejected before output', 'docs/alpha_phase_2026-10-27/PROOF.md', 'change'),
        ('changed old report rejected before output', 'experiments/alpha_phase_2026_10_27/exact_results.json', 'change'),
        ('missing imported kernel rejected before output', 'experiments/v04_w_interval/symkernel.py', 'remove'),
        ('extra unmanifested file rejected before output', 'unlisted.txt', 'extra')):
        fixture = output/('fixture-'+str(len(tests)))
        shutil.copytree(package, fixture)
        target = fixture/relative
        if action == 'remove':
            target.unlink()
        elif action == 'change':
            target.write_bytes(target.read_bytes()+b'\nINTENTIONAL REJECTION CONTROL\n')
        else:
            target.write_text('INTENTIONAL REJECTION CONTROL', encoding='utf-8')
        forbidden_output = output/('must-not-exist-'+str(len(tests)))
        def attempt():
            preflight(fixture, pin)
            new_output_directory(fixture, forbidden_output)
        refused(label, attempt)
        check(label+' leaves no output', not forbidden_output.exists())

    phase = package/'experiments/alpha_phase_2026_10_27'
    for name in REPORTS:
        original = read_json(phase/name)
        runtime = copy.deepcopy(original)
        runtime['elapsed_seconds'] += 1
        check('timing is non-scientific: '+name,
              canonical_digest(mathematical_payload(name, original)) == canonical_digest(mathematical_payload(name, runtime)))
        wrong = copy.deepcopy(original)
        wrong['status'] = 'INJECTED_WRONG_RESULT'
        check('mathematical/report result change detected: '+name,
              canonical_digest(mathematical_payload(name, original)) != canonical_digest(mathematical_payload(name, wrong)))

    bad_claims = read_json(package/'V050_CLAIMS.json')
    bad_claims['inner_inclusion_providers'] = ['v0.4.0']
    refused('misattributed universal inclusion rejected', lambda: validate_claims(bad_claims))
    malicious_zip = output/'path-escape.zip'
    with zipfile.ZipFile(malicious_zip, 'x') as archive:
        archive.writestr('../must-not-escape.txt', 'bad')
    refused('archive traversal rejected', lambda: safe_extract(malicious_zip, output/'malicious-extraction'))
    check('archive traversal made no external file', not (output/'must-not-escape.txt').exists())

    cli = package/'tools/reproduce_v050.py'
    candidate_output = output/'preflight-created-nothing'
    result = subprocess.run([sys.executable, '-B', '-S', str(cli), '--out-dir', str(candidate_output),
                              '--expected-manifest-sha256', pin, '--preflight-only'], capture_output=True, encoding='utf-8')
    check('CLI read-only preflight succeeds', result.returncode == 0 and not candidate_output.exists())
    existing_result = subprocess.run([sys.executable, '-B', '-S', str(cli), '--out-dir', str(existing)],
                                      capture_output=True, encoding='utf-8')
    check('CLI repeated output attempt fails before child execution', existing_result.returncode != 0 and
          digest(existing/'sentinel.txt') == before and not (existing/'raw').exists())
    verify_inventory(package, pin)
    write_new_json(output/'guard_tests.json', {'status': 'PASS', 'tests': tests, 'test_count': len(tests),
        'candidate_unchanged': True, 'manifest_sha256': pin,
        'qualification': 'Release-integrity controls only; not independent mathematical proof.'})
    print('PASS', len(tests), 'release-integrity controls; frozen candidate unchanged')


if __name__ == '__main__':
    main()
