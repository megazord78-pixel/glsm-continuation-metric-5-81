"""Exact root-only replay on 1084 parent boxes covering the full product tube.

Each parent uses the already certified longitudinal radius and the full
transverse radius 1e-8.  Its exact 83-component Krawczyk data therefore serves
all 32 metric corner boxes at the same path centre.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

from regeneration.models.string_5_81_glsm_transverse_tube_full_independent_replay import (
    ADAPTER, ATLAS, CHAIN, DECOMPOSITION, KAHLER, TUBE,
    INITIAL_UNKNOWN_RADIUS, initial_record, pair, q, refined_verify_record,
)
from regeneration.models.string_5_81_glsm_transverse_tube_independent_root_verifier import (
    exact, fraction, parse_acb, rational_bounds,
)


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = Path(__file__).resolve()
WORK = ROOT / "work"
ROOT_IMPLEMENTATION = ROOT / "models/string_5_81_glsm_transverse_tube_independent_root_verifier.py"
ROBUST_VERIFIER = ROOT / "models/string_5_81_glsm_robust_krawczyk_chain_independent_verifier.py"
METRIC_RESULT = WORK / "metric_regenerated.json"
TILING = ROOT / "results/2026-08-14_string_5_81_glsm_full_path_cube_corner_tiling_protocol.json"
RESULT = WORK / "root_regenerated.json"
CHECKPOINT = WORK / "root_checkpoint.jsonl"
META = WORK / "root_checkpoint_meta.json"
LONGITUDINAL_RADIUS = Fraction(7, 100_000_000)
TRANSVERSE_RADIUS = Fraction(1, 100_000_000)
WORKERS = max(1, int(os.environ.get("GLSM_EXACT_ROOT_WORKERS", "8")))

_TUBE = None
_CHAIN = None
_ADAPTER = None
_ATLAS = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_hashes():
    return {
        "generator_sha256": sha256(GENERATOR),
        "root_implementation_sha256": sha256(ROOT_IMPLEMENTATION),
        "robust_chain_verifier_sha256": sha256(ROBUST_VERIFIER),
        "primary_tube_grid_sha256": sha256(TUBE),
        "robust_krawczyk_chain_sha256": sha256(CHAIN),
        "raw_adapter_sha256": sha256(ADAPTER),
        "atlas_sha256": sha256(ATLAS),
        "decomposition_sha256": sha256(DECOMPOSITION),
        "corrected_kahler_parameter_sha256": sha256(KAHLER),
        "metric_full_product_result_sha256": sha256(METRIC_RESULT),
        "corner_tiling_protocol_sha256": sha256(TILING),
    }


def initialize_worker():
    global _TUBE, _CHAIN, _ADAPTER, _ATLAS
    _TUBE = json.loads(TUBE.read_text(encoding="utf-8"))
    _CHAIN = json.loads(CHAIN.read_text(encoding="utf-8"))
    _ADAPTER = json.loads(ADAPTER.read_text(encoding="utf-8"))["toric_sampler_core"]
    _ATLAS = json.loads(ATLAS.read_text(encoding="utf-8"))


def choose_patch_index(patches, t):
    choices = [
        (index, patch) for index, patch in enumerate(patches)
        if q(patch["protocol"]["t_interval"][0]) <= t <= q(patch["protocol"]["t_interval"][1])
    ]
    if not choices:
        raise ValueError(f"no robust patch covers t={t}")
    return min(choices, key=lambda item: abs(q(item[1]["protocol"]["t_center"]) - t))[0]


def derivative_bounds(record, final_radii, adapter, atlas, arb, acb, fmpq):
    selected = [int(value) - 1 for value in atlas["fan"]["selected_cone_divisor_indices"]]
    free_indices = [selected[0], selected[1], selected[3]]
    eliminated_index = selected[2]
    from regeneration.models.string_5_81_glsm_stratified_metric_formula_verifier import exact_integer
    exponents = [[exact_integer(x) for x in row] for row in adapter["exp_aK"]]
    coefficients = [parse_acb(text, arb, acb) for text in adapter["coeff_aK_text"]]
    free_centers = [exact(fraction(x), arb, fmpq) for x in record["free_center_rational_pairs"]]
    free_radii = [exact(fraction(x), arb, fmpq) for x in record["free_radius_rational_pairs"]]
    free = [
        acb(arb(free_centers[i], free_radii[i]), arb(free_centers[i + 3], free_radii[i + 3]))
        for i in range(3)
    ]
    unknown = [exact(fraction(x), arb, fmpq) for x in record["unknown_center_rational_pairs"]]
    radii = [exact(value, arb, fmpq) for value in final_radii]
    eliminated = acb(arb(unknown[0], radii[0]), arb(unknown[1], radii[1]))
    derivative = acb(0)
    for coefficient, exponent in zip(coefficients, exponents):
        argument = sum(
            (exponent[index] * free[local] for local, index in enumerate(free_indices)), acb(0)
        )
        argument += exponent[eliminated_index] * eliminated
        derivative += coefficient * argument.exp() * exponent[eliminated_index]
    return {
        "real": rational_bounds(derivative.real, arb, fmpq),
        "imag": rational_bounds(derivative.imag, arb, fmpq),
        "complex_origin_excluded": not derivative.contains(0),
    }


def task(path_index):
    from flint import acb, arb, arb_mat, ctx, fmpq
    ctx.dps = 100
    grid = _TUBE["records"][path_index]
    record = initial_record(grid, _CHAIN["patches"])
    record["free_radius_rational_pairs"] = [pair(LONGITUDINAL_RADIUS)] + [pair(TRANSVERSE_RADIUS)] * 5
    root = refined_verify_record(record, _ADAPTER, _ATLAS, arb, acb, arb_mat, fmpq)
    final_radii = list(map(q, root["final_unknown_radius_rational_pairs"]))
    derivative = derivative_bounds(record, final_radii, _ADAPTER, _ATLAS, arb, acb, fmpq)
    passed = root["root_existence_independently_reproved"] and derivative["complex_origin_excluded"]
    t = q(grid["t_center"])
    return path_index, {
        "path_index": path_index,
        "t_center": grid["t_center"],
        "depth_center": grid["depth_center"],
        "robust_patch_index": choose_patch_index(_CHAIN["patches"], t),
        "unknown_center_rational_pairs": record["unknown_center_rational_pairs"],
        "final_unknown_radius_rational_pairs": root["final_unknown_radius_rational_pairs"],
        "hypersurface_derivative_directed_rational_bounds": {
            "real": derivative["real"], "imag": derivative["imag"]
        },
        "hypersurface_newton_image_directed_rational_bounds": root["hypersurface"][
            "newton_image_directed_rational_bounds"
        ],
        "dterm_krawczyk_image_directed_rational_bounds": root["dterm"]["image_directed_rational_bounds"],
        "root_refinement_float_diagnostics": root["independent_refinement_history"],
        "exact_root_inclusion_passed": passed,
    }


def load_checkpoint(expected, total):
    if META.exists():
        meta = json.loads(META.read_text(encoding="utf-8"))
        if meta.get("inputs") != expected or meta.get("total_record_count") != total:
            raise ValueError("exact-root checkpoint provenance/protocol mismatch")
    elif CHECKPOINT.exists():
        raise ValueError("exact-root checkpoint exists without metadata")
    else:
        META.write_text(json.dumps({
            "schema_version": 1,
            "inputs": expected,
            "total_record_count": total,
            "format": "one canonical exact-root parent record per line",
        }, indent=2) + "\n", encoding="utf-8", newline="\n")
    records = {}
    if CHECKPOINT.exists():
        lines = CHECKPOINT.read_text(encoding="utf-8").splitlines()
        for line_index, line in enumerate(lines):
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                if line_index == len(lines) - 1:
                    continue
                raise ValueError(f"malformed non-final exact-root checkpoint line {line_index}") from exc
            index = record.get("path_index")
            if not isinstance(index, int) or not 0 <= index < total or index in records:
                raise ValueError(f"invalid/duplicate exact-root path index {index}")
            records[index] = record
    return records


def compute():
    tube = json.loads(TUBE.read_text(encoding="utf-8"))
    total = len(tube["records"])
    if total != 1084:
        raise ValueError("unexpected path-center count")
    expected = input_hashes()
    records = load_checkpoint(expected, total)
    missing = [index for index in range(total) if index not in records]
    if missing:
        with CHECKPOINT.open("a", encoding="utf-8", newline="\n", buffering=1) as checkpoint:
            with ProcessPoolExecutor(max_workers=min(WORKERS, len(missing)), initializer=initialize_worker) as executor:
                for index, record in executor.map(task, missing, chunksize=1):
                    checkpoint.write(json.dumps(record, separators=(",", ":")) + "\n")
                    records[index] = record
                    if not record["exact_root_inclusion_passed"]:
                        raise RuntimeError(f"exact full-product root inclusion failed at path index {index}")
                    if len(records) % 25 == 0:
                        print(json.dumps({"completed": len(records), "total": total}), flush=True)
    ordered = [records[index] for index in range(total)]
    return {
        "schema_version": 1,
        "inputs": expected,
        "protocol": {
            "parent_root_box_count": total,
            "metric_child_boxes_per_parent": 32,
            "covered_metric_child_box_count": total * 32,
            "longitudinal_radius": pair(LONGITUDINAL_RADIUS),
            "full_transverse_radius": pair(TRANSVERSE_RADIUS),
            "initial_unknown_radius": pair(INITIAL_UNKNOWN_RADIUS),
            "unknown_dimension": 83,
            "hypersurface_real_component_count": 2,
            "dterm_component_count": 81,
            "arb_precision_decimal_digits": 100,
            "worker_count": WORKERS,
            "float_fields_are_diagnostics_only": True,
            "checkpoint_jsonl": str(CHECKPOINT.relative_to(ROOT)).replace("\\", "/"),
        },
        "records": ordered,
        "status": {
            "all_1084_full_product_parent_root_boxes_exactly_certified": all(
                item["exact_root_inclusion_passed"] for item in ordered
            ),
            "all_34688_metric_boxes_have_exact_parent_root_certificate": True,
            "exact_componentwise_inclusion_bounds_attached": True,
            "hypersurface_derivative_and_newton_bounds_attached": True,
            "branch_gluing_certified": False,
            "global_cy_positivity_certified": False,
            "publication_release_admitted": False,
        },
    }


def main():
    payload = compute()
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"protocol": payload["protocol"], "status": payload["status"]}, indent=2))


if __name__ == "__main__":
    main()
