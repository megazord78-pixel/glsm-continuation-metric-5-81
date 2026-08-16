"""Stdlib-only exact verifier for the local six-real-dimensional open-set theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "data" / "endpoint_validation"
ARTIFACT = RESULTS / "2026-08-13_string_5_81_glsm_depth15_open_neighborhood_certificate.json"
FILES = {
    "nested_material_box_sha256": RESULTS / "2026-08-13_string_5_81_glsm_depth15_nested_material_interval_box.json",
    "independent_joint_root_sha256": RESULTS / "2026-08-13_string_5_81_glsm_depth15_independent_joint_root_krawczyk.json",
    "independent_baseline_metric_sha256": RESULTS / "2026-08-13_string_5_81_glsm_depth15_independent_baseline_metric.json",
    "independent_joint_correction_sha256": RESULTS / "2026-08-13_string_5_81_glsm_depth15_nested_independent_joint_verifier.json",
    "independent_materiality_consensus_sha256": RESULTS / "2026-08-13_string_5_81_glsm_depth15_independent_materiality_consensus.json",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def q(value) -> Fraction:
    if not isinstance(value, dict) or set(value) != {"numerator", "denominator"}:
        raise ValueError("noncanonical rational pair")
    if (not isinstance(value["numerator"], int) or isinstance(value["numerator"], bool)
            or not isinstance(value["denominator"], int) or isinstance(value["denominator"], bool)):
        raise ValueError("rational entries must be integers")
    if value["denominator"] <= 0:
        raise ValueError("nonpositive denominator")
    result = Fraction(value["numerator"], value["denominator"])
    if (result.numerator, result.denominator) != (value["numerator"], value["denominator"]):
        raise ValueError("unreduced rational pair")
    return result


def interval(entry, component):
    return q(entry[component]["lower"]), q(entry[component]["upper"])


def add(left, right):
    return left[0] + right[0], left[1] + right[1]


def neg(value):
    return -value[1], -value[0]


def mul(left, right):
    values = (left[0] * right[0], left[0] * right[1], left[1] * right[0], left[1] * right[1])
    return min(values), max(values)


def cmul(left, right):
    return add(mul(left[0], right[0]), neg(mul(left[1], right[1]))), add(
        mul(left[0], right[1]), mul(left[1], right[0])
    )


def load_matrix(serialized):
    return [[(interval(entry, "real"), interval(entry, "imag")) for entry in row] for row in serialized]


def determinant(matrix, dimension):
    zero = ((Fraction(0), Fraction(0)), (Fraction(0), Fraction(0)))
    one = ((Fraction(1), Fraction(1)), (Fraction(0), Fraction(0)))
    total = zero
    for permutation in itertools.permutations(range(dimension)):
        term = one
        for row, column in enumerate(permutation):
            term = cmul(term, matrix[row][column])
        inversions = sum(permutation[i] > permutation[j] for i in range(dimension) for j in range(i + 1, dimension))
        if inversions % 2:
            term = (neg(term[0]), neg(term[1]))
        total = (add(total[0], term[0]), add(total[1], term[1]))
    return total


def verify_payload(payload, check_live_provenance=True):
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported schema")
    if check_live_provenance:
        expected = {name: sha256(path) for name, path in FILES.items()}
        if any(payload["inputs"].get(name) != digest for name, digest in expected.items()):
            raise ValueError("live data provenance mismatch")
    nested = json.loads(FILES["nested_material_box_sha256"].read_text(encoding="utf-8"))
    root = json.loads(FILES["independent_joint_root_sha256"].read_text(encoding="utf-8"))
    baseline = json.loads(FILES["independent_baseline_metric_sha256"].read_text(encoding="utf-8"))
    correction = json.loads(FILES["independent_joint_correction_sha256"].read_text(encoding="utf-8"))
    materiality = json.loads(FILES["independent_materiality_consensus_sha256"].read_text(encoding="utf-8"))
    domain = payload["theorem_domain"]
    radius = q(domain["closed_box_common_free_radius"])
    unknown_radius = q(domain["implicit_unknown_radius"])
    center = list(map(q, domain["closed_box_center_rational_pairs"]))
    if domain["chart_complex_dimension"] != 3 or domain["chart_real_dimension"] != 6 or len(center) != 6:
        raise ValueError("chart dimension mismatch")
    if center != [Fraction(-3, 2)] + [Fraction(0)] * 5 or radius != Fraction(1, 10**9):
        raise ValueError("unexpected center or radius")
    if unknown_radius != Fraction(1, 10**8) or domain["implicit_unknown_real_dimension"] != 83:
        raise ValueError("unexpected implicit box")
    if list(map(q, nested["exact_center"]["free_center_rational_pairs"])) != center:
        raise ValueError("center mismatch")
    if {nested["protocol"]["free_radius"], root["protocol"]["free_radius"], correction["protocol"]["root_box_radii"]["free"]} != {"1e-9"}:
        raise ValueError("free radii mismatch")
    if {nested["protocol"]["eliminated_radius"], nested["protocol"]["dterm_radius"], root["protocol"]["eliminated_radius"], root["protocol"]["dterm_radius"], correction["protocol"]["root_box_radii"]["eliminated"], correction["protocol"]["root_box_radii"]["dterm"]} != {"1e-8"}:
        raise ValueError("implicit radii mismatch")
    records = root["component_records"]
    if len(records) != 83:
        raise ValueError("incomplete root inclusion")
    root_ratios = []
    for index, record in enumerate(records):
        upper, target = q(record["image_displacement_abs_upper_rational"]), q(record["target_radius_rational"])
        if record["index"] != index or not 0 <= upper < target or not record["strictly_inside"]:
            raise ValueError("root inclusion failed")
        root_ratios.append(upper / target)
    base_matrix = load_matrix(baseline["baseline_metric_entry_directed_rational_bounds"])
    correction_matrix = load_matrix(correction["correction_metric_entry_directed_rational_bounds"])
    corrected = [[(add(base_matrix[i][j][0], correction_matrix[i][j][0]), add(base_matrix[i][j][1], correction_matrix[i][j][1])) for j in range(3)] for i in range(3)]
    minor_lowers = []
    serialized = materiality["independent_corrected_metric_fraction_interval_minors"]
    for dimension, item in zip((1, 2, 3), serialized):
        real, imaginary = determinant(corrected, dimension)
        expected = (q(item["real_lower"]), q(item["real_upper"]), q(item["imaginary_lower"]), q(item["imaginary_upper"]))
        if item["dimension"] != dimension or expected != (*real, *imaginary) or real[0] <= 0 or not item["strictly_positive"]:
            raise ValueError("corrected Sylvester minor failed")
        minor_lowers.append(real[0])
    largest_materiality = max(q(record["absolute_relative_interval"][0]) for record in materiality["records"])
    if largest_materiality != q(payload["replayed_evidence"]["largest_relative_correction_lower"]):
        raise ValueError("materiality lower bound mismatch")
    if largest_materiality < Fraction(1, 100):
        raise ValueError("materiality threshold failed")
    depth = payload["path_embedding"]
    depth_center = q(depth["depth_center"])
    closed_path = depth["closed_box_one_sided_depth_interval"]
    open_path = depth["open_neighborhood_path_subsegment"]
    lower = q(closed_path["lower"])
    upper = q(closed_path["upper"])
    if depth_center != Fraction(3, 2) or (lower, upper) != (depth_center - radius, depth_center):
        raise ValueError("path subsegment is not exactly embedded")
    if (q(open_path["lower"]), q(open_path["upper"])) != (lower, upper):
        raise ValueError("open path bounds mismatch")
    if open_path["lower_endpoint_included"] is not False or open_path["upper_endpoint_included"] is not True:
        raise ValueError("open path endpoint topology mismatch")
    status = payload["status"]
    required_true = (
        "common_rational_box_exactly_matched_across_inputs",
        "unique_parameterized_joint_root_over_closed_box_certified",
        "corrected_metric_positive_on_closed_box_certified",
        "nonempty_open_six_real_dimensional_chart_neighborhood_certified",
        "corrected_metric_positive_on_open_neighborhood_certified",
        "open_neighborhood_contains_half_open_one_sided_path_subsegment",
        "independent_correction_materiality_at_least_one_percent_certified",
    )
    if not all(status.get(key) is True for key in required_true):
        raise ValueError("a local theorem status was dropped")
    for key in ("open_tube_along_entire_1_49985_to_1_5_path_certified", "global_corrected_metric_positivity_certified", "global_ricci_or_monge_ampere_error_bound_certified"):
        if status.get(key) is not False:
            raise ValueError("local theorem was promoted beyond evidence")
    if payload["replayed_evidence"] != {
        "joint_root_component_count": 83,
        "baseline_minor_count": 3,
        "corrected_minor_count": 3,
        "materiality_direction_count": len(materiality["records"]),
        "materiality_witness_count": materiality["summary"]["certified_witness_count"],
        "largest_relative_correction_lower": materiality["summary"]["largest_independent_relative_lower"],
    }:
        raise ValueError("evidence ledger mismatch")
    return {
        "free_radius_exact": str(radius),
        "open_real_dimension": 6,
        "joint_root_component_count": 83,
        "maximum_root_image_ratio_exact": str(max(root_ratios)),
        "corrected_metric_minor_lower_bounds_exact": [str(value) for value in minor_lowers],
        "largest_relative_correction_lower_exact": str(largest_materiality),
        "one_sided_path_depth_interval_exact": [str(lower), str(upper)],
        "nonempty_open_neighborhood_verified": True,
        "stdlib_exact_replay_passed": True,
    }


def main():
    print(json.dumps(verify_payload(json.loads(ARTIFACT.read_text(encoding="utf-8"))), indent=2))


if __name__ == "__main__":
    main()
