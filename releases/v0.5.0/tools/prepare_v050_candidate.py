"""Build/freeze a selected, unpublished v0.5.0 candidate without Git mutations.

All research/reference bytes remain original. This builder never runs the
unsafe historical aggregate. --repack works from an already frozen bundle.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

from v050_common import (REPORTS, NEW_SCOPE, OLD_SCOPE, canonical_digest, check_evidence_set,
    digest, mathematical_payload, read_json, require, verify_inventory, write_new_json)

ROOT = Path(__file__).resolve().parents[1]
AUDIT_DIR = 'scratch/alpha_independent_audit_2026-10-28'
PHASE_DIR = 'experiments/alpha_phase_2026_10_27'
PUBLISHED_SHA = '20328551989ac1fe7ff7bd5d5b6c28a9eac52bf5d7b8dcbd11b3ac92866481f9'


def claims():
    return {'candidate_id': 'v0.5.0-rc.1', 'release_status': 'UNPUBLISHED_RELEASE_CANDIDATE',
        'date': '2026-10-28', 'timezone': 'Asia/Shanghai', 'model_frozen': True,
        'universal_claim': 'intersection_{K>=6} W_K = W_6 intersect W_10 = [a10,b6]',
        'lower_binding_size': 10, 'upper_binding_size': 6,
        'inner_inclusion_providers': ['finite_cover', 'uniform_probe'],
        'baseline_alpha': ['1999/8000', '2001/8000'], 'exact_results_scope': NEW_SCOPE,
        'historical_scope_correction': {'json_pointer': '/scope', 'old': OLD_SCOPE, 'new': NEW_SCOPE},
        'stages': {
            'exact_results.json': {'stage': 'verify_phase', 'scope': NEW_SCOPE,
                'universal_inner_inclusion': False, 'shared_dependencies': ['canonical model', 'exact polynomial libraries', 'Theorem A']},
            'algebra_audit.json': {'stage': 'audit_phase',
                'scope': 'Formal moment identities, uniform prefix constants, and finite cross-checks; not sample-based universal positivity.',
                'shared_dependencies': ['some moment formulas', 'exact vector/polynomial primitives']},
            'finite_cover.json': {'stage': 'finite_cover',
                'scope': 'All real alpha in B=[123/500,7/24], all integer t>=0, K=7..203 except 10; plus K=10 on [1/4,7/24].',
                'infinite_time_basis': 'Analytic residence and quadratic tail certificate, with finite binary-word over-approximation.'},
            'finite_audit.json': {'stage': 'audit_cover',
                'scope': 'Reconstruct all 1046 finite cells using Sturm signs and four rejection controls.',
                'independence_limit': 'Second sign algorithm; shares model/trajectory algebra and analytic residence with the producer.'},
            'uniform_boxes.json': {'stage': 'uniform_probe',
                'scope': 'All K>=204, all real alpha in B=[123/500,7/24], all integer t>=0.',
                'infinite_domains': ['h=1/K in [0,1/204]', 'closed alpha/rounding boxes', 'analytic entire infinite tail'],
                'shared_dependencies': ['Theorem A', 'exact moment identities', 'rational interval arithmetic']}
        },
        'not_claimed': ['complete fixed-K phase diagram', 'general finite membership horizon',
                        'new mathematical research', 'human peer review', 'proof-assistant formalization', 'novelty'],
        'qualification': 'Reproduction and consolidation of the already-audited computer-assisted theorem; execution is not independent proof.'}


def create_zip(package, archive):
    require(not archive.exists(), 'archive exists; use a new path')
    archive.parent.mkdir(parents=True, exist_ok=True)
    manifest, pin = verify_inventory(package)
    names = sorted([*manifest['files'], 'V050_MANIFEST.json'])
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as bundle:
        for name in names:
            info = zipfile.ZipInfo('sand-v0.5.0-rc.1/' + name, date_time=(2026, 10, 28, 0, 0, 0))
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            bundle.writestr(info, (package/name).read_bytes(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return manifest, pin


def build(package):
    require(not package.exists(), 'candidate directory exists; choose a new one')
    audited = read_json(ROOT/AUDIT_DIR/'original_hashes.json')
    require(len(audited) == 699, 'unexpected original audit inventory')
    for name, expected in audited.items():
        require(digest(ROOT/name) == expected, 'audited workspace input changed: ' + name)
    private_manifest = read_json(ROOT/'docs/V040_MANIFEST.json')
    private_entries = {n: d for group in private_manifest['files'].values() for n, d in group.items()}
    for name, expected in private_entries.items():
        require(digest(ROOT/name) == expected, 'private v0.4.0 freeze changed: ' + name)
    reference_dir = ROOT/PHASE_DIR
    expected_payloads = {name: canonical_digest(mathematical_payload(name, read_json(reference_dir/name))) for name in REPORTS}
    check_evidence_set(ROOT, reference_dir, expected_payloads)
    for name, expected in read_json(reference_dir/'verification.json')['reports_sha256'].items():
        require(digest(reference_dir/name) == expected, 'original aggregate/report binding')

    package.mkdir(parents=True, exist_ok=False)
    tracked = {n for n in subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT, text=True).split('\0') if n}
    modified = set(subprocess.check_output(['git', 'diff', '--name-only'], cwd=ROOT, text=True).splitlines())
    git_head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    inventory = {}

    def put_bytes(name, data, role, origin):
        require(name not in inventory, 'duplicate selected destination: ' + name)
        target = package/name
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as stream:
            stream.write(data)
        inventory[name] = {'sha256': digest(target), 'bytes': len(data), 'role': role, 'origin': origin}

    def copy(source, role, destination=None):
        p = ROOT/source
        origin = {'type': 'workspace_file', 'source_path': source, 'source_sha256': digest(p),
                  'index_state': ('tracked_modified' if source in modified else 'tracked') if source in tracked else 'untracked',
                  'pinned_in_oct28_original_inventory': source in audited}
        if source in audited:
            require(origin['source_sha256'] == audited[source], 'audited bytes differ at copy time')
        put_bytes(destination or source, p.read_bytes(), role, origin)

    for directory, role in [
        (PHASE_DIR, 'original_alpha_sources_evidence_or_exploration'),
        ('docs/alpha_phase_2026-10-27', 'original_alpha_proof_or_figure'),
        ('docs/alpha_phase_audit_2026-10-28', 'original_oct28_audit'),
        ('docs/review_2026-10-27', 'preserved_review'),
        ('docs/corrections_2026-10-27', 'preserved_corrections')]:
        for p in sorted((ROOT/directory).rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts:
                copy(p.relative_to(ROOT).as_posix(), role)

    aggregate = read_json(reference_dir/'verification.json')
    required = set(aggregate['source_sha256']) | set(aggregate['prior_review_and_corrections_unchanged'])
    for name in REPORTS:
        required.update(read_json(reference_dir/name)['source_sha256'])
    required.update(('tools/audit_v04_interval_v3.py', 'tools/audit_v04_uniform_tail.py',
                     'tools/audit_v040_review.py', 'experiments/v04_w_interval/symkernel.py',
                     'docs/XO_ALPHA_EXTENSION_PROOF.md', 'docs/S3_XO_BOUNDARY_THEOREM_PROOF.md'))
    for name in sorted(required):
        if name not in inventory:
            copy(name, 'unchanged_math_dependency_or_original_source_binding')

    for name in ('independent_checks.py', 'independent_checks.json', 'audit_integrity.json',
                 'original_hashes.json', 'reproduction.log'):
        copy(AUDIT_DIR+'/'+name, 'original_oct28_audit_evidence')
    for name in (*REPORTS, 'verification.json', 'baseline_verification.json'):
        copy(AUDIT_DIR+'/reproduction/'+PHASE_DIR+'/'+name, 'original_oct28_isolated_reproduction_evidence')

    for name in sorted(private_entries):
        copy(name, 'private_research_v040_frozen_reference', 'references/research-v0.4.0/'+name)
    copy('docs/V040_MANIFEST.json', 'private_research_v040_frozen_manifest',
         'references/research-v0.4.0/docs/V040_MANIFEST.json')

    archive_source = 'release_records/v0.4.0/v0.4.0-verification-source.zip'
    require(digest(ROOT/archive_source) == PUBLISHED_SHA, 'published archive differs from preserved receipt')
    archive_dest = 'references/published-v0.4.0/v0.4.0-verification-source.zip'
    copy(archive_source, 'published_v040_archive', archive_dest)
    copy('release_records/v0.4.0/PUBLISHED.json', 'published_v040_receipt', 'references/published-v0.4.0/PUBLISHED.json')
    published = read_json(ROOT/'release_records/v0.4.0/PUBLISHED.json')
    asset = next(a for a in published['assets'] if a['name'] == 'v0.4.0-verification-source.zip')
    require(asset['digest'] == 'sha256:'+PUBLISHED_SHA, 'publication receipt/archive mismatch')
    with zipfile.ZipFile(ROOT/archive_source) as bundle:
        for name in ('LICENSE', 'LICENSE-CC-BY-4.0.txt', 'CONTENT_LICENSE.md'):
            member = 'v0.4.0-verification-source/'+name
            data = bundle.read(member)
            put_bytes(name, data, 'license_from_published_baseline', {'type': 'published_zip_member',
                'archive_sha256': PUBLISHED_SHA, 'member': member, 'source_sha256': hashlib.sha256(data).hexdigest()})

    for name in ('v050_common.py', 'reproduce_v050.py', 'prepare_v050_candidate.py', 'test_v050_reproduction.py'):
        copy('tools/'+name, 'new_release_engineering')
    copy('docs/v050_release/README.md', 'new_release_guide', 'README.md')
    copy('docs/v050_release/PROVENANCE.md', 'new_release_provenance', 'docs/V050_PROVENANCE.md')
    put_bytes('.gitattributes', b'* -text\n', 'new_release_byte_preservation',
              {'type': 'generated', 'purpose': 'Prevent future Git EOL conversion of hash-pinned bundle files.'})
    claim_blob = (json.dumps(claims(), ensure_ascii=False, indent=2)+'\n').encode('utf-8')
    put_bytes('V050_CLAIMS.json', claim_blob, 'new_release_scope_ledger',
              {'type': 'generated_metadata', 'generator': 'tools/prepare_v050_candidate.py', 'scientific_claim_changed': False})
    manifest = {'candidate_id': 'v0.5.0-rc.1', 'status': 'UNPUBLISHED_RELEASE_CANDIDATE',
        'prepared_date': '2026-10-28', 'timezone': 'Asia/Shanghai',
        'workspace_origin': {'git_head': git_head, 'head_is_not_complete_content_identity': True,
                             'untracked_research_materials_included_by_hash': True},
        'independent_audit': {'report': 'docs/alpha_phase_audit_2026-10-28/AUDIT_CN.md',
                              'original_inventory': AUDIT_DIR+'/original_hashes.json', 'all_699_originals_matched_before_build': True},
        'published_baseline': {'archive': archive_dest, 'sha256': PUBLISHED_SHA,
                               'root_directory': 'v0.4.0-verification-source', 'frozen_file_count': 78,
                               'scope': 'Only the published v0.4.0 narrow-interval theorem.'},
        'audited_payload_sha256': expected_payloads,
        'semantic_comparison_exceptions': {'all_five_reports': ['/elapsed_seconds'],
            'finite_audit.json': ['/witness_sha256, after checking the actual witness binding'],
            'consolidated_reports_only': ['/proof_scope addition', '/scope correction in exact_results.json; raw output retained']},
        'self_hash_policy': 'Manifest is the sole unlisted file; external receipt pins its SHA-256 and the archive.',
        'omissions': ['Git metadata and unrelated private research files',
                      'Large duplicate source tree from the October 28 isolation run; its linked reports are included',
                      'New reproduction outputs, which must be written to a separate new directory'],
        'files': dict(sorted(inventory.items()))}
    write_new_json(package/'V050_MANIFEST.json', manifest)
    # Exercise the same full preflight used by a clean extracted package.
    from reproduce_v050 import preflight
    preflight(package)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate-dir', type=Path, required=True)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--repack', action='store_true', help='Repack an already verified manifest; no original workspace or Git required')
    args = parser.parse_args()
    require(not args.receipt.exists() and not args.archive.exists(), 'archive/receipt already exists')
    package = args.candidate_dir.resolve()
    if not args.repack:
        build(package)
    manifest, pin = create_zip(package, args.archive.resolve())
    write_new_json(args.receipt, {'candidate_id': manifest['candidate_id'], 'status': 'FROZEN_NOT_PUBLISHED',
        'manifest_sha256': pin, 'archive_name': args.archive.name, 'archive_sha256': digest(args.archive),
        'archive_bytes': args.archive.stat().st_size, 'package_file_count_including_manifest': len(manifest['files'])+1,
        'deterministic_zip': {'ordered_members': True, 'timestamp': '2026-10-28T00:00:00', 'compression': 'DEFLATE level 9'},
        'scientific_sources_changed': False, 'original_reports_rewritten': False,
        'approval_to_tag_or_publish': False})
    print(json.dumps({'candidate_id': manifest['candidate_id'], 'files': len(manifest['files'])+1,
                      'manifest_sha256': pin, 'archive_sha256': digest(args.archive)}, indent=2))


if __name__ == '__main__':
    main()
