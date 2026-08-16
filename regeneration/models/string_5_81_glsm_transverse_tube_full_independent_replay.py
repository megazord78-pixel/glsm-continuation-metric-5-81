"""Parallel resumable raw-root + independent-metric replay of all tube boxes."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

from regeneration.models.string_5_81_glsm_transverse_tube_independent_root_verifier import (
    ADAPTER,
    ATLAS,
    DECOMPOSITION,
    KAHLER,
    pair,
    refined_verify_record,
)
from regeneration.models.string_5_81_glsm_transverse_tube_independent_formula_verifier import (
    PARAMETERS,
    verify_box,
)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
GENERATOR = Path(__file__).resolve()
ROOT_IMPLEMENTATION = ROOT / "models/string_5_81_glsm_transverse_tube_independent_root_verifier.py"
FORMULA_IMPLEMENTATION = ROOT / "models/string_5_81_glsm_transverse_tube_independent_formula_verifier.py"
TUBE = RESULTS / "2026-08-13_string_5_81_glsm_transverse_tube_cover.json"
CHAIN = RESULTS / "2026-08-13_string_5_81_glsm_robust_krawczyk_chain.json"
RESULT = RESULTS / "2026-08-13_string_5_81_glsm_transverse_tube_full_independent_replay.json"
CHECKPOINT = RESULTS / "2026-08-13_string_5_81_glsm_transverse_tube_full_independent_replay_checkpoint.json"
FREE_RADII = [Fraction("7e-8")] + [Fraction("1e-9")] * 5
INITIAL_UNKNOWN_RADIUS = Fraction("7e-7")
WORKERS = max(1, int(os.environ.get("GLSM_INDEPENDENT_TUBE_WORKERS", "6")))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def q(value):
    return Fraction(value["numerator"], value["denominator"])


def patch_contains(patch, t):
    lower, upper = map(q, patch["protocol"]["t_interval"])
    return lower <= t <= upper


def choose_patch(patches, t):
    choices = [patch for patch in patches if patch_contains(patch, t)]
    if not choices:
        raise ValueError(f"no affine predictor covers t={t}")
    return min(choices, key=lambda patch: abs(q(patch["protocol"]["t_center"]) - t))


def initial_record(tube_record, patches):
    t = q(tube_record["t_center"])
    patch = choose_patch(patches, t)
    local = patch["local_affine_coordinates"]
    slopes = list(map(q, local["slope_rational_pairs"]))
    intercepts = list(map(q, local["intercept_rational_pairs"]))
    unknown = [intercept + slope * t for intercept, slope in zip(intercepts, slopes)]
    free = [Fraction(-3, 2) - t] + [Fraction(0)] * 5
    return {
        "tube_record_index": tube_record["index"],
        "t_center": tube_record["t_center"],
        "depth_center": tube_record["depth_center"],
        "free_center_rational_pairs": [pair(value) for value in free],
        "free_radius_rational_pairs": [pair(value) for value in FREE_RADII],
        "unknown_center_rational_pairs": [pair(value) for value in unknown],
        "unknown_radius_rational_pairs": [pair(INITIAL_UNKNOWN_RADIUS)] * 83,
        "primary_leading_principal_minor_rational_lower_bounds": tube_record[
            "leading_principal_minor_rational_lower_bounds"
        ],
    }


def task(index):
    from flint import acb, acb_mat, arb, arb_mat, ctx, fmpq
    ctx.dps = 100
    tube = json.loads(TUBE.read_text(encoding="utf-8"))
    chain = json.loads(CHAIN.read_text(encoding="utf-8"))
    adapter = json.loads(ADAPTER.read_text(encoding="utf-8"))["toric_sampler_core"]
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    record = initial_record(tube["records"][index], chain["patches"])
    root = refined_verify_record(record, adapter, atlas, arb, acb, arb_mat, fmpq)
    if not root["root_existence_independently_reproved"]:
        return index, {
            "index": index,
            "independent_root_certified": False,
            "independent_metric_positive": False,
            "full_independent_replay_passed": False,
            "failure_stage": "independent root refinement",
            "root_refinement_history": root["independent_refinement_history"],
        }
    metric_record = dict(record)
    metric_record["unknown_radius_rational_pairs"] = root["final_unknown_radius_rational_pairs"]
    metric_record["primary_hypersurface_newton_image_directed_rational_bounds"] = root[
        "hypersurface"
    ]["newton_image_directed_rational_bounds"]
    metric_record["primary_dterm_krawczyk_image_directed_rational_bounds"] = root[
        "dterm"
    ]["image_directed_rational_bounds"]
    metric = verify_box(metric_record, arb, acb, acb_mat, arb_mat, fmpq)
    independent_minors = [
        item["directed_rational_bounds"] for item in metric["independent_leading_principal_minors"]
    ]
    passed = metric["all_three_independent_minors_positive"]
    return index, {
        "index": index,
        "t_center": tube["records"][index]["t_center"],
        "depth_center": tube["records"][index]["depth_center"],
        "root_refinement_history": root["independent_refinement_history"],
        "independent_root_dterm_maximum_image_ratio": root["dterm"]["maximum_image_radius_ratio"],
        "independent_joint_neumann_rho_upper": metric["joint_neumann_rho_upper"],
        "independent_moment_neumann_rho_upper": metric["moment_neumann_rho_upper"],
        "independent_leading_principal_minor_directed_rational_bounds": independent_minors,
        "independent_root_certified": True,
        "independent_metric_positive": passed,
        "full_independent_replay_passed": passed,
    }


def input_hashes():
    return {
        "generator_sha256": sha256(GENERATOR),
        "independent_root_implementation_sha256": sha256(ROOT_IMPLEMENTATION),
        "independent_formula_implementation_sha256": sha256(FORMULA_IMPLEMENTATION),
        "primary_tube_grid_sha256": sha256(TUBE),
        "affine_predictor_chain_sha256": sha256(CHAIN),
        "raw_adapter_sha256": sha256(ADAPTER),
        "atlas_sha256": sha256(ATLAS),
        "decomposition_sha256": sha256(DECOMPOSITION),
        "corrected_kahler_parameter_sha256": sha256(KAHLER),
        "exact_correction_parameters_sha256": sha256(PARAMETERS),
    }


def reusable(inputs):
    if not CHECKPOINT.exists():
        return {}
    payload = json.loads(CHECKPOINT.read_text(encoding="utf-8"))
    if payload.get("inputs") != inputs:
        return {}
    return {item["index"]: item for item in payload.get("records", [])}


def write_checkpoint(inputs, records):
    CHECKPOINT.write_text(json.dumps({
        "inputs": inputs,
        "records": [records[index] for index in sorted(records)],
    }, indent=2) + "\n", encoding="utf-8", newline="\n")


def compute():
    tube = json.loads(TUBE.read_text(encoding="utf-8"))
    inputs = input_hashes()
    records = reusable(inputs)
    missing = [index for index in range(len(tube["records"])) if index not in records]
    if missing:
        with ProcessPoolExecutor(max_workers=min(WORKERS, len(missing))) as executor:
            for index, record in executor.map(task, missing, chunksize=1):
                records[index] = record
                if len(records) % 20 == 0 or not record["full_independent_replay_passed"]:
                    write_checkpoint(inputs, records)
    ordered = [records[index] for index in range(len(tube["records"]))]
    passed = all(item["full_independent_replay_passed"] for item in ordered)
    minimum_third = min(
        q(item["independent_leading_principal_minor_directed_rational_bounds"][2]["lower"])
        for item in ordered if item["independent_metric_positive"]
    )
    return {
        "schema_version": 1,
        "inputs": inputs,
        "protocol": {
            "center_count": len(ordered),
            "uniform_initial_unknown_radius": pair(INITIAL_UNKNOWN_RADIUS),
            "free_radius_rational_pairs": [pair(value) for value in FREE_RADII],
            "independent_dterm_refinement_steps": 3,
            "independent_dterm_refinement_safety_factor": pair(Fraction(6, 5)),
            "raw_root_then_independent_metric_per_center": True,
            "primary_component_radii_not_reused": True,
            "worker_count": WORKERS,
            "resumable_checkpoint": str(CHECKPOINT.relative_to(ROOT)).replace("\\", "/"),
        },
        "records": ordered,
        "summary": {
            "independent_root_pass_count": sum(item["independent_root_certified"] for item in ordered),
            "independent_metric_pass_count": sum(item["independent_metric_positive"] for item in ordered),
            "minimum_independent_third_minor_lower_exact": str(minimum_third),
            "maximum_independent_root_q": max(item["independent_root_dterm_maximum_image_ratio"] for item in ordered),
        },
        "status": {
            "all_1084_root_boxes_independently_reproved": passed,
            "all_1084_metrics_formula_independently_reimplemented": passed,
            "complete_path_open_tube_independently_replayed": passed,
            "global_corrected_metric_positivity_on_cy_certified": False,
            "publication_release_admitted": False,
        },
    }


def refresh_provenance_only():
    """Update transitive implementation hashes after exact record replay."""
    from regeneration.models.string_5_81_glsm_transverse_tube_full_independent_replay_verifier import verify_payload

    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    verify_payload(payload, check_live_provenance=False)
    checkpoint = json.loads(CHECKPOINT.read_text(encoding="utf-8"))
    checkpoint_records = checkpoint.get("records", [])
    result_records = payload.get("records", [])
    if checkpoint_records != result_records[:len(checkpoint_records)]:
        raise ValueError("checkpoint is not an exact result prefix; provenance-only migration refused")
    fresh = input_hashes()
    payload["inputs"] = fresh
    checkpoint["inputs"] = fresh
    checkpoint["records"] = result_records
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    CHECKPOINT.write_text(json.dumps(checkpoint, indent=2) + "\n", encoding="utf-8", newline="\n")
    return {"verified_record_count": len(payload["records"]), "updated_inputs": fresh}


def main():
    payload = compute()
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"protocol": payload["protocol"], "summary": payload["summary"], "status": payload["status"]}, indent=2))


if __name__ == "__main__":
    main()
