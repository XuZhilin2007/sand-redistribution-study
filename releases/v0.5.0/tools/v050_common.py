"""Release-integrity utilities. No Sand mathematics is changed here."""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path, PurePosixPath
import stat
import zipfile

REPORTS = ('exact_results.json', 'algebra_audit.json', 'finite_cover.json',
           'finite_audit.json', 'uniform_boxes.json')
OLD_SCOPE = 'Only the exact theorems explicitly listed; universal inner interval inherited from v0.4.0'
NEW_SCOPE = ('This stage proves the listed low-dimensional fixed-K results and universal exclusions only. '
             'The new universal inner inclusion is supplied by finite_cover and uniform_probe, '
             'with the documented analytic dependencies; it is not inherited from v0.4.0.')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_new_json(path, data):
    path = Path(path)
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
        stream.write('\n')


def canonical_digest(data):
    blob = json.dumps(data, sort_keys=True, ensure_ascii=False, separators=(',', ':'), allow_nan=False)
    return hashlib.sha256(blob.encode('utf-8')).hexdigest()


def relative_name(name):
    require(isinstance(name, str) and name and '\\' not in name and ':' not in name, 'unsafe manifest name')
    p = PurePosixPath(name)
    require(not p.is_absolute() and all(x not in ('', '.', '..') for x in p.parts), 'unsafe relative path: ' + name)
    require(p.as_posix() == name, 'noncanonical path: ' + name)
    return p


def verify_inventory(root, expected_manifest_sha256=None):
    """Check the PRE-EXISTING freeze, including its exact complete file list."""
    root = Path(root).resolve()
    manifest_file = root/'V050_MANIFEST.json'
    actual_manifest = digest(manifest_file)
    if expected_manifest_sha256 is not None:
        require(actual_manifest == expected_manifest_sha256, 'manifest digest differs from external pin')
    manifest = read_json(manifest_file)
    require(manifest['candidate_id'] == 'v0.5.0-rc.1' and manifest['status'] == 'UNPUBLISHED_RELEASE_CANDIDATE',
            'unexpected candidate identity/status')
    entries = manifest['files']
    require(isinstance(entries, dict) and entries, 'empty inventory')
    expected = set(entries) | {'V050_MANIFEST.json'}
    require(len({n.casefold() for n in expected}) == len(expected), 'case-colliding inventory names')
    actual = set()
    for p in root.rglob('*'):
        require(not p.is_symlink(), 'symlink in frozen candidate: ' + str(p))
        if p.is_file():
            actual.add(p.relative_to(root).as_posix())
    require(actual == expected, 'candidate file list differs: missing=' + str(sorted(expected-actual)) +
            ', extra=' + str(sorted(actual-expected)))
    for name, item in entries.items():
        p = root/relative_name(name)
        require(p.resolve().is_relative_to(root), 'inventory escaped root')
        require(p.stat().st_size == item['bytes'] and digest(p) == item['sha256'], 'frozen file changed: ' + name)
    require(manifest['self_hash_policy'] == 'Manifest is the sole unlisted file; external receipt pins its SHA-256 and the archive.',
            'unexpected self-hash policy')
    return manifest, actual_manifest


def new_output_directory(root, out_dir):
    root, out = Path(root).resolve(), Path(out_dir).resolve()
    require(not out.is_relative_to(root), 'output directory must be outside the frozen candidate')
    require(not out.exists(), 'output directory already exists; refusing to overwrite evidence')
    out.mkdir(parents=True, exist_ok=False)
    return out


def mathematical_payload(name, record):
    """Only predeclared runtime fields are omitted, never mathematical fields."""
    require(name in REPORTS, 'unknown mathematical report')
    copied = dict(record)
    require(isinstance(copied.get('elapsed_seconds'), (int, float)), 'missing runtime field')
    del copied['elapsed_seconds']
    if name == 'finite_audit.json':
        value = copied.pop('witness_sha256')
        require(isinstance(value, str) and len(value) == 64, 'invalid witness binding')
    return copied


def verify_source_bindings(root, reports):
    root = Path(root).resolve()
    for report in reports.values():
        for name, expected in report['source_sha256'].items():
            p = root/relative_name(name.replace('\\', '/'))
            require(p.resolve().is_relative_to(root) and digest(p) == expected, 'source binding failed: ' + name)


def check_evidence_set(root, directory, expected_payloads):
    directory = Path(directory)
    reports = {name: read_json(directory/name) for name in REPORTS}
    verify_source_bindings(root, reports)
    require(reports['exact_results.json']['scope'] == OLD_SCOPE, 'unexpected historical scope; investigate')
    require(reports['finite_audit.json']['witness_sha256'] == digest(directory/'finite_cover.json'),
            'finite audit is not bound to the actual witness')
    for name, record in reports.items():
        require(canonical_digest(mathematical_payload(name, record)) == expected_payloads[name],
                'mathematical output differs from audited evidence: ' + name)
    return reports


def assemble(reports):
    """Retain the audited aggregate's explicit root/order/coverage checks."""
    exact, algebra, finite, second, boxes = (reports[n] for n in REPORTS)
    roots = exact['central']['roots']
    a6, aa, bb = [[Fraction(x) for x in roots[name]['bracket']] for name in ('a6', 'a10', 'b6')]
    require(a6[1] < Fraction(123, 500) < aa[0] < aa[1] < Fraction(1, 4) < bb[0] < bb[1] < Fraction(7, 24),
            'exact endpoint order')
    require(roots['a10']['irreducible_mod_prime'] == 29 and roots['b6']['irreducible_mod_prime'] == 11,
            'boundary degrees')
    targets = set(range(7, 204))-{10}
    require(set(finite['targets']) == set(second['verified_targets']) == targets and
            len(finite['targets']) == len(second['verified_targets']) == len(targets), 'incomplete finite target set')
    require(second['K10_upper_verified'], 'missing K10 upper coverage')
    require(algebra['formal_identity_count'] == 60 and len(algebra['exact_vector_cases']) == 95, 'algebra coverage')
    require(boxes['status'] == 'ALL_BOXES_POSITIVE' and boxes['failure_count'] == 0, 'uniform proof failed')
    require(boxes['K_min'] == 204 and boxes['alpha'] == ['123/500', '7/24'], 'uniform domain')
    require(boxes['branch_alpha'] == {'M': ['123/500', '13/50'], 'U': ['1/4', '7/24']}, 'branch domain')
    for name, floor in boxes['strict_lower_bounds'].items():
        require(Fraction(boxes['minima'][name]['lower']) > Fraction(floor), 'strict margin ' + name)
    return {'boundary_roots': {name: roots[name] for name in ('a10', 'b6')},
            'finite_K_count': len(targets), 'finite_cells_including_K10': second['verified_cells'],
            'formal_identity_count': algebra['formal_identity_count'],
            'box_count': boxes['box_count'], 'terminal_box_count': boxes['terminal_box_count']}


def safe_extract(archive, destination):
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    with zipfile.ZipFile(archive) as bundle:
        seen = set()
        for entry in bundle.infolist():
            name = entry.filename.rstrip('/')
            p = relative_name(name)
            require(name.casefold() not in seen, 'duplicate archive path')
            seen.add(name.casefold())
            require(not stat.S_ISLNK(entry.external_attr >> 16), 'symlink archive member')
            target = destination/p
            require(target.resolve().is_relative_to(destination), 'archive path escape')
            if entry.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open('xb') as stream:
                    stream.write(bundle.read(entry))

