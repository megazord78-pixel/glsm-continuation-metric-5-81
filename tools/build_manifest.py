"""Create the deterministic immutable release manifest."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "MANIFEST.sha256.json"
EXCLUDED_DIRS = {"__pycache__", ".pytest_cache", "build", "dist", ".git", "work", "tmp"}
EXCLUDED_SUFFIXES = {".pyc", ".pyo", ".aux", ".blg", ".log", ".out", ".toc"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    files = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file()
        and path != OUTPUT
        and not any(part in EXCLUDED_DIRS for part in path.relative_to(ROOT).parts)
        and path.suffix not in EXCLUDED_SUFFIXES
    )
    payload = {
        "schema_version": 1,
        "release_version": "2.0.0",
        "hash_algorithm": "SHA-256",
        "excluded": {
            "directories": sorted(EXCLUDED_DIRS),
            "suffixes": sorted(EXCLUDED_SUFFIXES),
            "reason": "ephemeral caches, regeneration checkpoints, PDF-QA renders, and compiler intermediates",
        },
        "files": {
            path.relative_to(ROOT).as_posix(): sha256(path)
            for path in files
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"manifest_created": True, "file_count": len(files)}, indent=2))


if __name__ == "__main__":
    main()
