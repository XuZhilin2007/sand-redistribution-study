"""Minimal test-assertion counter for the sand redistribution study.

Run:  python tests/count_tests.py

Executes every tests/test_*.py, counts runtime [PASS] lines (the project's
assertion-count convention), and prints per-suite and total counts. Added
2026-09-19 after three completion-report counting errors; use this instead
of manual arithmetic in completion reports.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    tests_dir = Path(__file__).resolve().parent
    files = sorted(tests_dir.glob("test_*.py"))
    total = 0
    failures = []
    for f in files:
        proc = subprocess.run([sys.executable, str(f)], capture_output=True, text=True)
        passed = proc.stdout.count("[PASS]")
        failed = proc.stdout.count("[FAIL]")
        total += passed
        status = "OK" if (failed == 0 and proc.returncode == 0) else "FAILING"
        if status != "OK":
            failures.append(f.name)
        print(f"{f.name}: {passed} assertions ({status})")
    print(f"TOTAL runtime assertions: {total}")
    if failures:
        print(f"FAILING suites: {failures}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
