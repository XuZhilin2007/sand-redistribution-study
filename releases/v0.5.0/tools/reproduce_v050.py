"""Safe v0.5.0 RC reproduction. Outputs must go to a NEW external directory.

Original algorithms and reference evidence are never rewritten. Exact
mathematical output is compared with PRE-FROZEN audited reference payloads.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

from v050_common import (REPORTS, OLD_SCOPE, NEW_SCOPE, assemble, canonical_digest,
    check_evidence_set, digest, new_output_directory, read_json, require,
    safe_extract, verify_inventory, verify_source_bindings, write_new_json)

ROOT = Path(__file__).resolve().parents[1]


def validate_claims(claims):
    require(claims['model_frozen'] is True and claims['release_status'] == 'UNPUBLISHED_RELEASE_CANDIDATE', 'claim identity')
    require(claims['universal_claim'] == 'intersection_{K>=6} W_K = W_6 intersect W_10 = [a10,b6]', 'claim changed')
    require(claims['exact_results_scope'] == NEW_SCOPE, 'incorrect exact-stage scope')
    require(claims['inner_inclusion_providers'] == ['finite_cover', 'uniform_probe'], 'incorrect inclusion provenance')
    require(claims['baseline_alpha'] == ['1999/8000', '2001/8000'], 'baseline scope expanded')
    require(set(claims['stages']) == set(REPORTS), 'incomplete stage scope ledger')


def preflight(root, pin=None):
    manifest, manifest_hash = verify_inventory(root, pin)
    claims = read_json(root/'V050_CLAIMS.json')
    validate_claims(claims)
    reference = root/'experiments/alpha_phase_2026_10_27'
    reports = check_evidence_set(root, reference, manifest['audited_payload_sha256'])
    original_aggregate = read_json(reference/'verification.json')
    verify_source_bindings(root, {'aggregate': original_aggregate})
    for name, expected in original_aggregate['reports_sha256'].items():
        require(digest(reference/name) == expected, 'original report was changed before reproduction: ' + name)
    for name, expected in original_aggregate['prior_review_and_corrections_unchanged'].items():
        require(digest(root/name) == expected, 'preserved review/correction binding failed: ' + name)
    require(digest(root/manifest['published_baseline']['archive']) == manifest['published_baseline']['sha256'], 'published baseline asset changed')
    return manifest, manifest_hash, claims, assemble(reports)


def run_stage(label, arguments, cwd, out):
    """Child stages get explicit output paths; -S and empty PATH test isolation."""
    command = [sys.executable, '-B', '-S', *map(str, arguments)]
    environment = {**os.environ, 'PATH': '', 'PYTHONDONTWRITEBYTECODE': '1', 'PYTHONNOUSERSITE': '1',
                   'PYTHONUTF8': '1'}
    environment.pop('PYTHONPATH', None)
    environment.pop('PYTHONHOME', None)
    tick = time.monotonic()
    print('Checking ' + label, flush=True)
    log = out/'logs'/f'{label}.log'
    with log.open('x', encoding='utf-8', newline='\n') as stream:
        process = subprocess.Popen(command, cwd=cwd, env=environment,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding='utf-8')
        for line in process.stdout:
            stream.write(line)
            stream.flush()
            print('  ' + line.rstrip(), flush=True)
        returncode = process.wait()
    record = {'stage': label, 'command': command, 'exit_code': returncode,
              'seconds': time.monotonic()-tick, 'log_sha256': digest(log)}
    require(returncode == 0, 'stage failed; see ' + str(log))
    return record


def reproduce(root, output, pin=None):
    started = time.monotonic()
    # Everything here is checked BEFORE creating an output directory or
    # spawning a process. Newly generated hashes cannot replace this gate.
    manifest, manifest_hash, claims, expected_summary = preflight(root, pin)
    out = new_output_directory(root, output)
    try:
        for name in ('raw', 'reports', 'logs'):
            (out/name).mkdir()
        write_new_json(out/'preflight.json', {'status': 'PASS', 'manifest_sha256': manifest_hash,
            'candidate_file_count': len(manifest['files'])+1, 'original_evidence_prechecked': True,
            'external_manifest_pin_checked': pin is not None,
            'proof_documents_frozen': ['docs/alpha_phase_2026-10-27/PROOF.md', 'docs/alpha_phase_2026-10-27/RESULTS.md']})
        source = root/'experiments/alpha_phase_2026_10_27'
        raw = out/'raw'
        commands = [
            ('verify_phase', [source/'verify_phase.py', '--report', raw/'exact_results.json']),
            ('audit_phase', [source/'audit_phase.py', '--report', raw/'algebra_audit.json']),
            ('finite_cover', [source/'finite_cover.py', '--max-k', '203', '--report', raw/'finite_cover.json']),
            ('audit_cover', [source/'audit_cover.py', '--input', raw/'finite_cover.json', '--report', raw/'finite_audit.json']),
            ('uniform_probe', [source/'uniform_probe.py', '--alpha-parts', '32', '--h-parts', '4', '--K-min', '204',
                               '--report', raw/'uniform_boxes.json'])]
        stages = [run_stage(label, args, root, out) for label, args in commands]
        regenerated = check_evidence_set(root, raw, manifest['audited_payload_sha256'])
        summary = assemble(regenerated)
        require(summary == expected_summary, 'assembled results differ from audited reference')
        # Correct only this new metadata; retain the historical producer's
        # exact raw output and its old scope value as historical evidence.
        corrections = []
        for name in REPORTS:
            record = dict(regenerated[name])
            if name == 'exact_results.json':
                require(record['scope'] == OLD_SCOPE, 'unexpected scope discrepancy')
                record['scope'] = NEW_SCOPE
                corrections.append({'report': name, 'field': '/scope', 'old': OLD_SCOPE, 'new': NEW_SCOPE,
                                    'classification': 'proof-scope metadata only; mathematical payload unchanged'})
            record['proof_scope'] = claims['stages'][name]
            roundtrip = dict(record)
            del roundtrip['proof_scope']
            if name == 'exact_results.json':
                roundtrip['scope'] = OLD_SCOPE
            require(roundtrip == regenerated[name], 'metadata consolidation changed a scientific field')
            write_new_json(out/'reports'/name, record)
        write_new_json(out/'metadata_corrections.json', {'corrections': corrections,
            'raw_reports_retained': True, 'scientific_fields_changed': []})

        published = manifest['published_baseline']
        safe_extract(root/published['archive'], out/'baseline_extract')
        baseline_root = out/'baseline_extract'/published['root_directory']
        require(not (baseline_root/'.git').exists(), 'baseline unexpectedly has Git metadata')
        stages.append(run_stage('published_v040', [baseline_root/'tools/verify_v040.py', '--report', raw/'baseline_verification.json'],
                                baseline_root, out))
        baseline = read_json(raw/'baseline_verification.json')
        require(baseline['status'] == 'PASS' and len(baseline['results']) == 6 and
                baseline['frozen_inputs_unchanged'] and baseline['frozen_file_count'] == published['frozen_file_count'],
                'published baseline incomplete')

        # The already-audited independent reconstruction writes next to its
        # script; execute an unchanged COPY inside the new output tree.
        audit_source = root/'scratch/alpha_independent_audit_2026-10-28/independent_checks.py'
        audit_copy = out/'independent_checks.py'
        shutil.copyfile(audit_source, audit_copy)
        stages.append(run_stage('oct28_audit_replay', [audit_copy], out, out))
        audit_result = read_json(out/'independent_checks.json')
        require(audit_result == read_json(root/'scratch/alpha_independent_audit_2026-10-28/independent_checks.json'),
                'October 28 audit reconstruction changed')

        verify_inventory(root, manifest_hash)
        # Recheck the original references AFTER execution, independently of
        # the newly written output metadata.
        preflight(root, manifest_hash)
        output_files = {p.relative_to(out).as_posix(): digest(p) for p in out.rglob('*')
                        if p.is_file() and not p.is_relative_to(out/'baseline_extract')}
        result = {'candidate_id': manifest['candidate_id'], 'status': 'REPRODUCED_AUDITED_CANDIDATE',
            'manifest_sha256': manifest_hash, 'candidate_files_unchanged': True,
            'all_mathematical_payloads_match_audited_reference': True,
            'proof_scope_metadata_corrected_in_new_reports_only': True,
            'universal_claim': claims['universal_claim'], 'summary': summary,
            'published_baseline': {'archive_sha256': published['sha256'], 'status': baseline['status'],
                                   'checks': 6, 'frozen_files': baseline['frozen_file_count'],
                                   'scope': 'v0.4.0 original narrow interval only; not the new B inclusion'},
            'oct28_independent_audit_replay_matches': True, 'stages': stages,
            'output_sha256_before_aggregate': output_files, 'seconds': time.monotonic()-started,
            'environment': {'python': sys.version, 'platform': platform.platform(), 'child_PATH': '',
                            'site_packages_disabled': True, 'git_required': False, 'network_required': False},
            'qualification': ('This reproduces the audited computation and its documented analytic assembly. '
                              'Successful execution is not an independent proof, formalization, human peer review, or publication approval.')}
        write_new_json(out/'verification.json', result)
        print(result['status'] + '; frozen candidate unchanged; audited payloads match', flush=True)
        return result
    except Exception as error:
        failure = out/'FAILURE.json'
        if not failure.exists():
            write_new_json(failure, {'status': 'FAILED', 'error': str(error),
                                    'original_evidence_overwritten': False})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, required=True, help='NEW directory outside the frozen package; existing paths are rejected')
    parser.add_argument('--expected-manifest-sha256', help='external digest from the candidate receipt')
    parser.add_argument('--preflight-only', action='store_true', help='check the complete freeze and old evidence without creating output')
    args = parser.parse_args()
    if args.preflight_only:
        preflight(ROOT, args.expected_manifest_sha256)
        require(not args.out_dir.resolve().is_relative_to(ROOT.resolve()) and not args.out_dir.exists(),
                'unsafe or existing output directory')
        print('PREFLIGHT PASS; no files written', flush=True)
    else:
        reproduce(ROOT, args.out_dir, args.expected_manifest_sha256)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as error:
        print('REFUSED/FAILED: ' + str(error), file=sys.stderr)
        raise SystemExit(1)
