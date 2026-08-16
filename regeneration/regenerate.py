"""Fresh Arb reconstruction and comparison with the archived release.

The regeneration layer is intentionally separate from the dependency-free
last-mile verifier.  It recreates all 34,688 metric child enclosures and all
1,084 parent root enclosures from frozen formula inputs using python-flint.
Checkpoint files live under ``regeneration/work`` and are not proof inputs.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
WORK = Path(__file__).resolve().parent / "work"


def load(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def require_environment():
    try:
        import flint
        from flint import ctx
    except ImportError as exc:
        raise SystemExit("full regeneration requires python-flint==0.9.0") from exc
    version = getattr(flint, "__version__", "unknown")
    if version != "0.9.0":
        raise SystemExit(f"expected python-flint 0.9.0, found {version}")
    ctx.dps = 100
    return version


def root_projection(record):
    return {
        key: record[key] for key in (
            "path_index", "t_center", "depth_center", "robust_patch_index",
            "unknown_center_rational_pairs", "final_unknown_radius_rational_pairs",
            "hypersurface_derivative_directed_rational_bounds",
            "hypersurface_newton_image_directed_rational_bounds",
            "dterm_krawczyk_image_directed_rational_bounds",
            "exact_root_inclusion_passed",
        )
    }


def metric_projection(record):
    return {
        key: record[key] for key in (
            "sequence", "path_index", "offset_slots", "t_center", "depth_center",
            "root_pass", "minor_directed_rational_bounds", "metric_pass", "full_pass",
        )
    }


def compare_records(name, regenerated, archived, projection):
    if len(regenerated) != len(archived):
        raise SystemExit(f"{name} record-count mismatch")
    for index, (left, right) in enumerate(zip(regenerated, archived)):
        if projection(left) != projection(right):
            raise SystemExit(f"{name} proof-field mismatch at record {index}")


def regenerate_all(verify_against_release: bool):
    version = require_environment()
    WORK.mkdir(parents=True, exist_ok=True)

    # Imports occur after the environment and work directory are established;
    # worker processes can therefore import the same standalone module graph.
    from regeneration.models import string_5_81_glsm_full_path_transverse_cube_cached_replay as metric
    from regeneration.models import string_5_81_glsm_full_product_exact_root_replay as root

    metric_payload = metric.compute()
    metric.RESULT.write_text(json.dumps(metric_payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    root_payload = root.compute()
    root.RESULT.write_text(json.dumps(root_payload, indent=2) + "\n", encoding="utf-8", newline="\n")

    if verify_against_release:
        archived_metric = load(PACKAGE / "data" / "metric_certificate.json.gz")
        archived_root = load(PACKAGE / "data" / "root_certificate.json.gz")
        compare_records("metric", metric_payload["records"], archived_metric["records"], metric_projection)
        compare_records("root", root_payload["records"], archived_root["records"], root_projection)

    print(json.dumps({
        "full_arb_regeneration_passed": True,
        "python_flint_version": version,
        "metric_child_boxes": len(metric_payload["records"]),
        "parent_root_boxes": len(root_payload["records"]),
        "proof_fields_equal_to_release": verify_against_release,
        "diagnostic_float_fields_compared": False,
    }, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-against-release", action="store_true")
    args = parser.parse_args()
    regenerate_all(args.verify_against_release)


if __name__ == "__main__":
    main()
