"""Verify the frozen v0.4.0 scientific snapshot without rewriting its inputs.

Uses the Python standard library. Checks the fixed-file manifest, runs the
five mathematical auditors and the affected regressions, then checks that
all frozen bytes are unchanged. --report optionally saves this run separately.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "docs/V040_MANIFEST.json"
CHECKS = (
    ("finite_common", "tools/audit_v04_interval_v3.py", "--scope", "common"),
    ("prefix", "tools/audit_v04_prefix_lemma.py"),
    ("W0", "tools/audit_v04_w0_uniform.py"),
    ("D2", "tools/audit_v04_d2_entry.py"),
    ("tail", "tools/audit_v04_uniform_tail.py"),
    ("regressions", "experiments/v04_w_interval/negative_tests_w.py"),
)


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest_files():
    record = json.loads(MANIFEST.read_text(encoding="utf-8"))
    expected = {}
    for entries in record["files"].values():
        for name, digest in entries.items():
            path = (ROOT / name).resolve()
            if ROOT not in path.parents:
                raise ValueError("manifest path outside snapshot: " + name)
            if name in expected and expected[name] != digest:
                raise ValueError("conflicting manifest entries: " + name)
            expected[name] = digest
    if not expected:
        raise ValueError("empty snapshot manifest")
    return expected


def check_hashes(expected):
    failures = []
    for name, digest in expected.items():
        path = ROOT / name
        if not path.is_file():
            failures.append({"file": name, "problem": "missing"})
        elif hash_file(path) != digest:
            failures.append({"file": name, "problem": "SHA-256 mismatch"})
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, help="save this run outside the frozen input set")
    args = parser.parse_args()
    started = time.monotonic()
    expected = manifest_files()
    protected = {(ROOT / name).resolve() for name in expected} | {MANIFEST.resolve()}
    if args.report and args.report.resolve() in protected:
        parser.error("--report would overwrite a frozen input; choose a new run path")
    manifest_sha = hash_file(MANIFEST)
    failures = check_hashes(expected)
    results = []
    if not failures:
        print(f"Snapshot SHA-256: {len(expected)} files match", flush=True)
        for name, script, *extra in CHECKS:
            print("Checking " + name + "...", flush=True)
            tick = time.monotonic()
            result = subprocess.run([sys.executable, "-B", script, *extra], cwd=ROOT,
                                    env={**os.environ, "PYTHONUTF8": "1"},
                                    capture_output=True, text=True, encoding="utf-8")
            results.append({"name": name, "command": ["python", "-B", script, *extra],
                            "exit_code": result.returncode,
                            "elapsed_seconds": round(time.monotonic() - tick, 3),
                            "stdout": result.stdout, "stderr": result.stderr})
            lines = result.stdout.splitlines()
            print(f"  exit {result.returncode}: " + (lines[-1] if lines else result.stderr[-500:]),
                  flush=True)
            if result.returncode:
                break
    after = check_hashes(expected)
    unchanged = not after and hash_file(MANIFEST) == manifest_sha
    passed = not failures and len(results) == len(CHECKS) and all(
        r["exit_code"] == 0 for r in results) and unchanged
    record = {"version": "v0.4.0", "status": "PASS" if passed else "FAIL",
              "manifest_sha256": manifest_sha, "frozen_file_count": len(expected),
              "initial_hash_failures": failures, "final_hash_failures": after,
              "frozen_inputs_unchanged": unchanged, "results": results,
              "elapsed_seconds": round(time.monotonic() - started, 3)}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if failures:
        print(json.dumps(failures, ensure_ascii=False))
    print(f"v0.4.0 {record['status']}; {len(results)}/{len(CHECKS)} checks; "
          f"frozen inputs unchanged={unchanged}; {record['elapsed_seconds']:.2f}s", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
