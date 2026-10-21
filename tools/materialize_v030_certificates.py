"""Restore the exact v0.3 B/E JSON certificates from checked gzip parts.

Uses only the Python standard library. Run once in a downloaded public
snapshot before running the original, unchanged scientific verifiers.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGETS = {
    "experiments/v03_b_lock_allk/results/b_lock_cert.json",
    "experiments/v03_e_cycle_allk/results/e_cycle_cert.json",
}
BLOCK = 1024 * 1024


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        while data := stream.read(BLOCK):
            value.update(data)
    return value.hexdigest()


def main() -> int:
    manifest = json.loads((ROOT / "certificates/v0.3.0/manifest.json").read_text(
        encoding="utf-8"))
    science = json.loads((ROOT / "docs/V030_FILE_MANIFEST.json").read_text(
        encoding="utf-8"))
    expected = {entry["path"]: entry for entry in science["files"]}
    files = manifest["files"]
    if manifest.get("format") != "gzip-parts-v1" or len(files) != 2 or \
            {entry["target"] for entry in files} != TARGETS:
        raise ValueError("Unexpected certificate transport set or format")
    for entry in files:
        rel = entry["target"]
        if entry["uncompressed_sha256"] != expected[rel]["sha256"] or \
                entry["uncompressed_size"] != expected[rel]["size"]:
            raise ValueError(f"Transport SHA does not match scientific manifest: {rel}")
        if not 0 < entry["uncompressed_size"] <= 512 * 1024 * 1024:
            raise ValueError(f"Unexpected certificate size: {rel}")
        target = (ROOT / rel).resolve()
        if not target.is_relative_to(ROOT.resolve()):
            raise ValueError("Certificate target escapes checkout")
        if target.exists():
            if sha256(target) != expected[rel]["sha256"]:
                raise ValueError(f"Preserving different existing certificate: {rel}")
            print("VERIFIED existing", rel, flush=True)
            continue
        parts = entry["parts"]
        if not parts or len({p["path"] for p in parts}) != len(parts):
            raise ValueError("Missing or duplicate certificate parts")
        for part in parts:
            path = (ROOT / part["path"]).resolve()
            if not path.is_relative_to((ROOT / "certificates/v0.3.0").resolve()):
                raise ValueError("Certificate part escapes transport directory")
            if path.stat().st_size != part["size"] or sha256(path) != part["sha256"]:
                raise ValueError(f"Compressed part hash/size mismatch: {part['path']}")
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(
                    mode="wb", prefix=".v030-", suffix=".tmp",
                    dir=target.parent, delete=False) as output:
                temporary = Path(output.name)
                value = hashlib.sha256()
                size = 0
                for part in parts:
                    with gzip.open(ROOT / part["path"], "rb") as stream:
                        while data := stream.read(BLOCK):
                            size += len(data)
                            if size > entry["uncompressed_size"]:
                                raise ValueError("Decompressed size exceeds manifest")
                            value.update(data)
                            output.write(data)
            if size != entry["uncompressed_size"] or \
                    value.hexdigest() != expected[rel]["sha256"]:
                raise ValueError(f"Restored certificate hash/size mismatch: {rel}")
            temporary.replace(target)
            temporary = None
            print("RESTORED", rel, "SHA256", expected[rel]["sha256"], flush=True)
        finally:
            if temporary is not None:
                if not temporary.resolve().is_relative_to(ROOT.resolve()):
                    raise ValueError("Unexpected temporary file")
                temporary.unlink(missing_ok=True)
    print("PASS: original scientific certificate bytes restored and verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
