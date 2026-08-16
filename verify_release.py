from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "MANIFEST.sha256.json"
EXCLUDED_DIRS = {"__pycache__", ".pytest_cache", "build", "dist", ".git", "work", "tmp"}
EXCLUDED_SUFFIXES = {".pyc", ".pyo", ".aux", ".blg", ".log", ".out", ".toc"}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    expected = manifest["files"]
    actual = sorted(
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file()
        and path.name != MANIFEST.name
        and not any(part in EXCLUDED_DIRS for part in path.relative_to(ROOT).parts)
        and path.suffix not in EXCLUDED_SUFFIXES
    )
    if actual != sorted(expected):
        raise SystemExit("release file set differs from manifest")
    for relative, digest in expected.items():
        if sha256(ROOT / relative) != digest:
            raise SystemExit(f"SHA-256 mismatch: {relative}")
    return len(expected)


def run(command, env=None):
    subprocess.run(command, cwd=ROOT, env=env, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="run manifest, exact last-mile verifier and tests")
    parser.add_argument("--full", action="store_true", help="also regenerate all Arb root and metric enclosures")
    args = parser.parse_args()
    full = args.full
    count = verify_manifest()
    env = os.environ.copy()
    source = str(ROOT / "src")
    env["PYTHONPATH"] = source + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    run([sys.executable, "-m", "glsmtube"], env=env)
    run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"], env=env)
    if full:
        run([sys.executable, "-m", "regeneration.regenerate", "--verify-against-release"], env=env)
    print(json.dumps({
        "release_manifest_verified": True,
        "quick_exact_verifier_passed": True,
        "packaged_mutation_tests_passed": True,
        "full_arb_regeneration_passed": full,
        "file_count": count,
    }, indent=2))


if __name__ == "__main__":
    main()
