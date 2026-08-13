from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "MANIFEST.sha256.json"
EXCLUDED_DIRS = {"__pycache__", ".pytest_cache", "build", "dist", ".git"}
EXCLUDED_SUFFIXES = {".pyc", ".pyo", ".aux", ".blg", ".log", ".out", ".toc"}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    expected = manifest["files"]
    actual_paths = sorted(
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*") if path.is_file()
        and not any(part in EXCLUDED_DIRS for part in path.parts)
        and path.name != MANIFEST.name
        and path.suffix not in EXCLUDED_SUFFIXES
    )
    if actual_paths != sorted(expected):
        raise SystemExit("release file set differs from manifest")
    for relative, digest in expected.items():
        if sha256(ROOT / relative) != digest:
            raise SystemExit(f"SHA-256 mismatch: {relative}")
    subprocess.run([sys.executable, "-m", "models.string_5_81_glsm_continuation_chain_independent_verifier"],
                   cwd=ROOT, check=True)
    subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"],
                   cwd=ROOT, check=True)
    print(json.dumps({"release_manifest_verified": True,
                      "central_verifier_passed": True,
                      "packaged_tests_passed": True,
                      "file_count": len(expected)}, indent=2))


if __name__ == "__main__":
    main()
