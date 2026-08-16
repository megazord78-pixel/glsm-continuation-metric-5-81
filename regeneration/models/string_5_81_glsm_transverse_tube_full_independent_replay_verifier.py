"""Stdlib-only exact verifier of the full independent tube replay."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
ARTIFACT = RESULTS / "2026-08-13_string_5_81_glsm_transverse_tube_full_independent_replay.json"
TUBE = RESULTS / "2026-08-13_string_5_81_glsm_transverse_tube_cover.json"
CHAIN = RESULTS / "2026-08-13_string_5_81_glsm_robust_krawczyk_chain.json"
GENERATOR = ROOT / "models/string_5_81_glsm_transverse_tube_full_independent_replay.py"
ROOT_IMPLEMENTATION = ROOT / "models/string_5_81_glsm_transverse_tube_independent_root_verifier.py"
FORMULA_IMPLEMENTATION = ROOT / "models/string_5_81_glsm_transverse_tube_independent_formula_verifier.py"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def q(value):
    if not isinstance(value, dict) or set(value) != {"numerator", "denominator"}:
        raise ValueError("noncanonical rational pair")
    result = Fraction(value["numerator"], value["denominator"])
    if value["denominator"] <= 0 or result.denominator != value["denominator"]:
        raise ValueError("noncanonical rational pair")
    return result


def verify_payload(payload, check_live_provenance=True):
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported schema")
    if check_live_provenance:
        if payload["inputs"]["generator_sha256"] != sha256(GENERATOR):
            raise ValueError("generator provenance mismatch")
        if payload["inputs"]["independent_root_implementation_sha256"] != sha256(ROOT_IMPLEMENTATION):
            raise ValueError("independent root implementation provenance mismatch")
        if payload["inputs"]["independent_formula_implementation_sha256"] != sha256(FORMULA_IMPLEMENTATION):
            raise ValueError("independent formula implementation provenance mismatch")
        if payload["inputs"]["primary_tube_grid_sha256"] != sha256(TUBE):
            raise ValueError("tube grid provenance mismatch")
        if payload["inputs"]["affine_predictor_chain_sha256"] != sha256(CHAIN):
            raise ValueError("predictor provenance mismatch")
    protocol = payload["protocol"]
    if protocol["center_count"] != 1084 or q(protocol["uniform_initial_unknown_radius"]) != Fraction(7, 10_000_000):
        raise ValueError("protocol geometry changed")
    if protocol["independent_dterm_refinement_steps"] != 3:
        raise ValueError("refinement depth changed")
    if q(protocol["independent_dterm_refinement_safety_factor"]) != Fraction(6, 5):
        raise ValueError("refinement safety factor changed")
    if not protocol["raw_root_then_independent_metric_per_center"] or not protocol["primary_component_radii_not_reused"]:
        raise ValueError("independence protocol dropped")
    tube = json.loads(TUBE.read_text(encoding="utf-8"))
    records = payload["records"]
    if len(records) != 1084 or len(tube["records"]) != 1084:
        raise ValueError("record count changed")
    minimum_third = None
    maximum_q = 0.0
    for index, (record, primary) in enumerate(zip(records, tube["records"])):
        if record["index"] != index or record["t_center"] != primary["t_center"] or record["depth_center"] != primary["depth_center"]:
            raise ValueError(f"grid mismatch at {index}")
        history = record["root_refinement_history"]
        if len(history) != 4 or [item["step"] for item in history] != [0, 1, 2, 3]:
            raise ValueError(f"refinement history mismatch at {index}")
        if any(item["included_component_count"] != 81 or item["maximum_image_radius_ratio"] >= 1 for item in history):
            raise ValueError(f"independent root inclusion failed at {index}")
        minors = record["independent_leading_principal_minor_directed_rational_bounds"]
        if len(minors) != 3:
            raise ValueError(f"minor dimension failed at {index}")
        for bounds in minors:
            lower, upper = q(bounds["lower"]), q(bounds["upper"])
            if not 0 < lower <= upper:
                raise ValueError(f"independent positivity failed at {index}")
        if not record["independent_root_certified"] or not record["independent_metric_positive"] or not record["full_independent_replay_passed"]:
            raise ValueError(f"record status failed at {index}")
        third = q(minors[2]["lower"])
        minimum_third = third if minimum_third is None else min(minimum_third, third)
        maximum_q = max(maximum_q, record["independent_root_dterm_maximum_image_ratio"])
    summary = payload["summary"]
    if summary["independent_root_pass_count"] != 1084 or summary["independent_metric_pass_count"] != 1084:
        raise ValueError("summary counts failed")
    if Fraction(summary["minimum_independent_third_minor_lower_exact"]) != minimum_third:
        raise ValueError("minimum third minor summary mismatch")
    if summary["maximum_independent_root_q"] != maximum_q:
        raise ValueError("maximum q summary mismatch")
    status = payload["status"]
    if not all(status[key] for key in (
        "all_1084_root_boxes_independently_reproved",
        "all_1084_metrics_formula_independently_reimplemented",
        "complete_path_open_tube_independently_replayed",
    )):
        raise ValueError("independent tube conclusion dropped")
    if status["global_corrected_metric_positivity_on_cy_certified"] or status["publication_release_admitted"]:
        raise ValueError("scope overclaimed")
    return {
        "verified_record_count": 1084,
        "verified_root_refinement_step_count": 4336,
        "verified_rational_determinant_interval_count": 3252,
        "minimum_independent_third_minor_lower_exact": str(minimum_third),
        "maximum_independent_root_q": maximum_q,
        "complete_path_open_tube_full_independent_replay_verified": True,
        "global_corrected_metric_positivity_on_cy": False,
        "stdlib_exact_replay_passed": True,
    }


def main():
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    print(json.dumps(verify_payload(payload), indent=2))


if __name__ == "__main__":
    main()
