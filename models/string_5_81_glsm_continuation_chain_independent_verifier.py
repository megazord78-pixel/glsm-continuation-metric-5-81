"""Independent stdlib-only verifier for the continuation theorem.

No continuation, metric, bridge, Arb, NumPy, or project generator is imported.
The verifier checks exact rational proof records and live byte hashes only.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/2026-08-13_string_5_81_glsm_continuation_chain_publication_certificate.json"

COVERS = (
    "2026-08-13_string_5_81_glsm_affine_path_metric_connected_cover.json",
    "2026-08-13_string_5_81_glsm_recentered_patch_metric_cover.json",
    "2026-08-13_string_5_81_glsm_recentered_patch_2_metric_cover.json",
    "2026-08-13_string_5_81_glsm_local_tangent_metric_cover.json",
    "2026-08-13_string_5_81_glsm_local_tangent_2_metric_cover.json",
    "2026-08-13_string_5_81_glsm_local_tangent_3_metric_cover.json",
)
BRIDGES = (
    "2026-08-13_string_5_81_glsm_root_patch_branch_gluing.json",
    "2026-08-13_string_5_81_glsm_root_patch_branch_gluing_2.json",
    "2026-08-13_string_5_81_glsm_local_tangent_branch_gluing.json",
    "2026-08-13_string_5_81_glsm_local_tangent_branch_gluing_2.json",
)
BRIDGE_3_COMPONENTS = (
    "2026-08-13_string_5_81_glsm_local_tangent_bridge_3_component_certificate.json"
)
ROBUST_CHAIN = "2026-08-13_string_5_81_glsm_robust_krawczyk_chain.json"
BASE_ROOT = "2026-08-13_string_5_81_glsm_affine_path_krawczyk_ladder.json"
MATERIAL_BOX = "2026-08-13_string_5_81_glsm_depth15_material_interval_box.json"
NESTED_MATERIAL_BOX = "2026-08-13_string_5_81_glsm_depth15_nested_material_interval_box.json"
TRANSITIVE_PROVENANCE_AUDIT = "2026-08-13_string_5_81_glsm_continuation_transitive_provenance_audit.json"
INDEPENDENT_MATERIAL_ROOT = "2026-08-13_string_5_81_glsm_depth15_independent_joint_root_krawczyk.json"
INDEPENDENT_BASELINE = "2026-08-13_string_5_81_glsm_depth15_independent_baseline_metric.json"
INDEPENDENT_CORRECTION = "2026-08-13_string_5_81_glsm_depth15_nested_independent_joint_verifier.json"
INDEPENDENT_MATERIALITY = "2026-08-13_string_5_81_glsm_depth15_independent_materiality_consensus.json"
STRATIFIED_METRIC_FORMULA = "2026-08-13_string_5_81_glsm_stratified_metric_formula_verifier.json"
EXPECTED_INTERVAL = (Fraction(29997, 20000), Fraction(3, 2))
EXPECTED_T_INTERVAL = (Fraction(-3, 20000), Fraction(0))
ROBUST_TARGET = Fraction(17, 20)
MATERIALITY_TARGET = Fraction(1, 100)
PUBLICATION_PACKAGE_NAME = "glsm-continuation-metric-5-81-v1.0.0"
GAUSSIAN_PHASES = ((1, 0), (-1, 0), (0, 1), (0, -1))


def publication_bundle_structurally_ready():
    # In the parent workspace the package is below publication/. In the
    # isolated release, ROOT itself is the package root.
    candidates = (ROOT / "publication" / PUBLICATION_PACKAGE_NAME, ROOT)
    required = (
        "LICENSE", "CITATION.cff", "README.md", "verify_release.py",
        "MANIFEST.sha256.json", "manuscript/continuation_theorem.md",
    )
    for package in candidates:
        if all((package / name).is_file() for name in required):
            return package, True
    return candidates[0], False


def fraction(pair: dict) -> Fraction:
    if set(pair) != {"numerator", "denominator"}:
        raise ValueError("noncanonical rational pair")
    value = Fraction(int(pair["numerator"]), int(pair["denominator"]))
    if value.denominator != int(pair["denominator"]):
        raise ValueError("unreduced rational pair")
    return value


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pair(value: Fraction) -> dict:
    return {"numerator": value.numerator, "denominator": value.denominator}


def _interval(bounds: dict) -> tuple[Fraction, Fraction]:
    lower, upper = fraction(bounds["lower"]), fraction(bounds["upper"])
    if lower > upper:
        raise ValueError("reversed rational interval")
    return lower, upper


def _add(a, b):
    return a[0] + b[0], a[1] + b[1]


def _neg(a):
    return -a[1], -a[0]


def _mul(a, b):
    values = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return min(values), max(values)


def _complex_add(a, b):
    return _add(a[0], b[0]), _add(a[1], b[1])


def _complex_mul(a, b):
    return (_add(_mul(a[0], b[0]), _neg(_mul(a[1], b[1]))),
            _add(_mul(a[0], b[1]), _mul(a[1], b[0])))


def _load_complex_matrix(serialized):
    if len(serialized) != 3 or any(len(row) != 3 for row in serialized):
        raise ValueError("metric matrix is not 3x3")
    return [[(_interval(entry["real"]), _interval(entry["imag"])) for entry in row]
            for row in serialized]


def _determinant(matrix, dimension):
    zero = ((Fraction(0), Fraction(0)), (Fraction(0), Fraction(0)))
    one = ((Fraction(1), Fraction(1)), (Fraction(0), Fraction(0)))
    total = zero
    for permutation in itertools.permutations(range(dimension)):
        term = one
        for row, column in enumerate(permutation):
            term = _complex_mul(term, matrix[row][column])
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(dimension) for j in range(i + 1, dimension))
        if inversions % 2:
            term = (_neg(term[0]), _neg(term[1]))
        total = _complex_add(total, term)
    return total


def _materiality_directions():
    result = []
    for index in range(3):
        vector = [(0, 0)] * 3
        vector[index] = (1, 0)
        result.append((f"axis_{index}", tuple(vector)))
    for left in range(3):
        for right in range(left + 1, 3):
            for phase_index, phase in enumerate(GAUSSIAN_PHASES):
                vector = [(0, 0)] * 3
                vector[left], vector[right] = (1, 0), phase
                result.append((f"pair_{left}_{right}_phase_{phase_index}", tuple(vector)))
    for second_index, second in enumerate(GAUSSIAN_PHASES):
        for third_index, third in enumerate(GAUSSIAN_PHASES):
            result.append((f"triple_phase_{second_index}_{third_index}", ((1, 0), second, third)))
    return tuple(result)


def _quadratic(matrix, vector):
    total = (Fraction(0), Fraction(0))
    for row in range(3):
        ar, ai = map(Fraction, vector[row])
        total = _add(total, tuple((ar * ar + ai * ai) * value for value in matrix[row][row][0]))
        for column in range(row + 1, 3):
            br, bi = map(Fraction, vector[column])
            real_coefficient = 2 * (ar * br + ai * bi)
            imag_coefficient = -2 * (ar * bi - ai * br)
            for coefficient, component in ((real_coefficient, matrix[row][column][0]),
                                           (imag_coefficient, matrix[row][column][1])):
                values = coefficient * component[0], coefficient * component[1]
                total = _add(total, (min(values), max(values)))
    return total


def verify_local_material_chain(root: dict, baseline: dict, correction: dict,
                                materiality: dict, hashes: dict,
                                endpoint_root: dict | None = None,
                                material_box: dict | None = None,
                                nested_box: dict | None = None) -> dict:
    if root.get("schema_version") != 1 or len(root.get("component_records", [])) != 83:
        raise ValueError("independent material root is incomplete")
    maximum_root_ratio = Fraction(0)
    for expected_index, record in enumerate(root["component_records"]):
        image = fraction(record["image_displacement_abs_upper_rational"])
        radius = fraction(record["target_radius_rational"])
        if record.get("index") != expected_index or image < 0 or not image < radius or not record.get("strictly_inside"):
            raise ValueError("independent material-root inclusion failed")
        maximum_root_ratio = max(maximum_root_ratio, image / radius)
    root_bounds = root.get("root_image_directed_rational_bounds", [])
    if len(root_bounds) != 83:
        raise ValueError("independent material root lacks 83 image intervals")
    for bounds in root_bounds:
        _interval(bounds)
    if not root["status"].get("independent_parameterized_joint_root_existence_and_uniqueness_certified"):
        raise ValueError("independent material-root theorem is false")
    if root["status"].get("primary_newton_or_krawczyk_images_reused"):
        raise ValueError("material root reuses primary images")

    endpoint_gluing = None
    if endpoint_root is not None or material_box is not None or nested_box is not None:
        if endpoint_root is None or material_box is None or nested_box is None:
            raise ValueError("endpoint gluing inputs must be supplied together")
        if root["inputs"].get("nested_center_sha256") != hashes[NESTED_MATERIAL_BOX]:
            raise ValueError("independent root is not bound to the live nested endpoint box")
        if material_box.get("exact_center") != nested_box.get("exact_center"):
            raise ValueError("material endpoint centers differ")
        if Fraction(str(material_box["protocol"].get("target_depth"))) != Fraction(3, 2):
            raise ValueError("material box is not at exact depth 3/2")
        if Fraction(str(nested_box["protocol"].get("target_depth"))) != Fraction(3, 2):
            raise ValueError("nested box is not at exact depth 3/2")
        if Fraction(str(endpoint_root["protocol"].get("center_depth"))) != Fraction(3, 2):
            raise ValueError("continuation endpoint root uses another depth")
        endpoint_record = next((record for record in endpoint_root.get("records", [])
                                if record.get("path_radius") == "1e-6"), None)
        if endpoint_record is None or not endpoint_record.get("certified"):
            raise ValueError("certified continuation endpoint record is absent")
        selected = endpoint_record.get("selected") or {}
        if not selected.get("full_parameterized_root_theorem"):
            raise ValueError("continuation endpoint existence theorem is absent")
        old_radii = selected.get("residual_component_radius_rational_pairs")
        if old_radii is None:
            scalar = selected.get("residual_component_radius_rational_pair")
            if scalar is None:
                raise ValueError("continuation endpoint radii are absent")
            old_radii = [scalar] * 83
        if len(old_radii) != 83:
            raise ValueError("continuation endpoint lacks 83 component radii")
        new_radii = [fraction(record["target_radius_rational"])
                     for record in root["component_records"]]
        old_radii_exact = list(map(fraction, old_radii))
        if any(radius <= 0 for radius in old_radii_exact):
            raise ValueError("continuation endpoint radius is nonpositive")
        if any(old > new for old, new in zip(old_radii_exact, new_radii)):
            raise ValueError("continuation endpoint root box is not nested in uniqueness box")
        endpoint_gluing = {
            "exact_depth": pair(Fraction(3, 2)),
            "componentwise_old_root_box_inside_new_uniqueness_box_count": 83,
            "maximum_old_to_new_radius_ratio_exact": pair(max(
                old / new for old, new in zip(old_radii_exact, new_radii)
            )),
            "same_exact_center_serialization": True,
            "old_endpoint_root_exists": True,
            "new_box_root_is_unique": True,
            "endpoint_root_identity_by_nested_existence_and_uniqueness": True,
        }

    if baseline["inputs"].get("independent_joint_root_sha256") != hashes[INDEPENDENT_MATERIAL_ROOT]:
        raise ValueError("baseline is not bound to the independent material root")
    if correction["inputs"].get("independent_joint_root_certificate_sha256") != hashes[INDEPENDENT_MATERIAL_ROOT]:
        raise ValueError("correction is not bound to the independent material root")
    if materiality["inputs"].get("independent_baseline_metric_sha256") != hashes[INDEPENDENT_BASELINE]:
        raise ValueError("materiality is not bound to the independent baseline")
    if materiality["inputs"].get("independent_joint_verifier_sha256") != hashes[INDEPENDENT_CORRECTION]:
        raise ValueError("materiality is not bound to the independent correction")

    baseline_matrix = _load_complex_matrix(baseline["baseline_metric_entry_directed_rational_bounds"])
    correction_matrix = _load_complex_matrix(correction["correction_metric_entry_directed_rational_bounds"])
    corrected = [[_complex_add(baseline_matrix[row][column], correction_matrix[row][column])
                  for column in range(3)] for row in range(3)]
    recomputed_minors = []
    serialized_minors = materiality.get("independent_corrected_metric_fraction_interval_minors", [])
    if len(serialized_minors) != 3:
        raise ValueError("local corrected metric lacks three serialized minors")
    for dimension, serialized in zip((1, 2, 3), serialized_minors):
        real, imaginary = _determinant(corrected, dimension)
        expected = (fraction(serialized["real_lower"]), fraction(serialized["real_upper"]),
                    fraction(serialized["imaginary_lower"]), fraction(serialized["imaginary_upper"]))
        if serialized.get("dimension") != dimension or expected != (*real, *imaginary):
            raise ValueError("local corrected minor does not recompute")
        if real[0] <= 0 or not serialized.get("strictly_positive"):
            raise ValueError("local corrected metric is not positive definite")
        recomputed_minors.append(real[0])

    directions = _materiality_directions()
    records = materiality.get("records", [])
    if len(records) != len(directions) or fraction(materiality["protocol"]["target_relative_lower_bound"]) != MATERIALITY_TARGET:
        raise ValueError("materiality protocol changed")
    witnesses = []
    largest = Fraction(0)
    for record, (label, vector) in zip(records, directions):
        serialized_vector = tuple((item["real"], item["imag"])
                                  for item in record["component_gaussian_integer_pairs"])
        if record.get("label") != label or serialized_vector != vector:
            raise ValueError("materiality direction rule changed")
        numerator = _quadratic(correction_matrix, vector)
        denominator = _quadratic(baseline_matrix, vector)
        if denominator[0] <= 0:
            raise ValueError("materiality denominator is nonpositive")
        absolute = (-numerator[1], -numerator[0]) if numerator[1] < 0 else (
            numerator if numerator[0] > 0 else (Fraction(0), max(-numerator[0], numerator[1])))
        ratio = absolute[0] / denominator[1], absolute[1] / denominator[0]
        if tuple(map(fraction, record["independent_correction_quadratic_interval"])) != numerator:
            raise ValueError("materiality correction quadratic mismatch")
        if tuple(map(fraction, record["independent_baseline_quadratic_interval"])) != denominator:
            raise ValueError("materiality baseline quadratic mismatch")
        if tuple(map(fraction, record["absolute_relative_interval"])) != ratio:
            raise ValueError("materiality relative interval mismatch")
        certified = (numerator[1] < 0 or numerator[0] > 0) and ratio[0] >= MATERIALITY_TARGET
        if record.get("relative_correction_at_least_one_percent_certified") != certified:
            raise ValueError("materiality witness status mismatch")
        largest = max(largest, ratio[0])
        if certified:
            witnesses.append(label)
    if materiality["summary"].get("certified_witness_labels") != witnesses:
        raise ValueError("materiality witness summary mismatch")
    status = materiality["status"]
    if not status.get("local_independent_corrected_metric_positive_definite_certified") or not witnesses:
        raise ValueError("local material theorem is not admitted")
    if status.get("global_corrected_metric_positivity_certified") or status.get("global_ricci_or_monge_ampere_error_bound_certified"):
        raise ValueError("local material theorem was promoted globally")
    result = {
        "joint_real_unknown_dimension": 83,
        "componentwise_root_inclusion_count": 83,
        "maximum_root_image_radius_ratio_exact": pair(maximum_root_ratio),
        "corrected_metric_minimum_leading_principal_minor_lower_bounds": [pair(value) for value in recomputed_minors],
        "materiality_direction_count": len(records),
        "materiality_witness_labels": witnesses,
        "largest_relative_correction_lower_bound_exact": pair(largest),
        "largest_relative_correction_lower_bound": float(largest),
        "local_corrected_metric_positive_definite": True,
        "global_metric_or_ricci_claim": False,
    }
    if endpoint_gluing is not None:
        result["continuation_endpoint_gluing"] = endpoint_gluing
    return result


def _segment_interval(segment: dict) -> tuple[Fraction, Fraction]:
    if "depth_lower" in segment:
        return fraction(segment["depth_lower"]), fraction(segment["depth_upper"])
    return fraction(segment["physical_depth_lower"]), fraction(segment["physical_depth_upper"])


def _strict_component_records(uppers, radii, label: str) -> Fraction:
    if len(uppers) != 83 or len(radii) != 83:
        raise ValueError(f"{label} lacks 83 rational component bounds")
    maximum = Fraction(0)
    for upper_pair, radius_pair in zip(uppers, radii):
        upper = fraction(upper_pair)
        radius = fraction(radius_pair)
        if upper < 0 or radius <= 0 or not upper < radius:
            raise ValueError(f"{label} component is not strictly included")
        maximum = max(maximum, upper / radius)
    return maximum


def verify_robust_chain(robust: dict, base_root_sha256: str,
                        material_box_sha256: str) -> dict:
    if robust.get("schema_version") != 1:
        raise ValueError("unsupported robust-chain schema")
    if robust["inputs"].get("base_root_sha256") != base_root_sha256:
        raise ValueError("robust chain is not anchored to the live base root")
    if robust["inputs"].get("material_box_sha256") != material_box_sha256:
        raise ValueError("robust chain is not anchored to the live exact center")
    if fraction(robust["protocol"]["target_maximum_image_radius_ratio"]) != ROBUST_TARGET:
        raise ValueError("robust contraction target changed")
    patches = robust.get("patches", [])
    bridges = robust.get("bridges", [])
    if len(patches) != 7 or len(bridges) != 7:
        raise ValueError("expected seven robust patches and seven bridges")
    intervals = []
    maximum_patch_ratio = Fraction(0)
    for patch in patches:
        lower, upper = map(fraction, patch["protocol"]["t_interval"])
        if not lower < upper:
            raise ValueError("degenerate robust patch")
        intervals.append((lower, upper))
        maximum_patch_ratio = max(maximum_patch_ratio, _strict_component_records(
            patch["krawczyk"]["component_image_abs_upper_rational_pairs"],
            patch["krawczyk"]["residual_component_radius_rational_pairs"],
            "robust patch",
        ))
    intervals.sort()
    if intervals[0][0] != EXPECTED_T_INTERVAL[0] or intervals[-1][1] != EXPECTED_T_INTERVAL[1]:
        raise ValueError("robust root-chain endpoints changed")
    if any(right[0] > left[1] for left, right in zip(intervals, intervals[1:])):
        raise ValueError("gap in robust root-chain cover")
    maximum_bridge_ratio = Fraction(0)
    for bridge in bridges:
        maximum_bridge_ratio = max(maximum_bridge_ratio, _strict_component_records(
            bridge["krawczyk"]["component_image_abs_upper_rational_pairs"],
            bridge["krawczyk"]["component_target_radius_rational_pairs"],
            "robust bridge",
        ))
    if maximum_patch_ratio >= ROBUST_TARGET or maximum_bridge_ratio >= ROBUST_TARGET:
        raise ValueError("robust chain does not meet its contraction target")
    return {
        "patch_count": len(patches),
        "bridge_count": len(bridges),
        "componentwise_patch_inclusion_count": len(patches) * 83,
        "componentwise_bridge_inclusion_count": len(bridges) * 83,
        "certified_t_interval": [pair(value) for value in EXPECTED_T_INTERVAL],
        "maximum_patch_image_radius_ratio_exact": pair(maximum_patch_ratio),
        "maximum_patch_image_radius_ratio": float(maximum_patch_ratio),
        "maximum_bridge_image_radius_ratio_exact": pair(maximum_bridge_ratio),
        "maximum_bridge_image_radius_ratio": float(maximum_bridge_ratio),
        "target_ratio_strictly_met": True,
    }


def verify_documents(covers: list[dict], bridges: list[dict], bridge3: dict,
                     robust: dict | None = None, base_root_sha256: str | None = None,
                     material_box_sha256: str | None = None) -> dict:
    if len(covers) != 6 or len(bridges) != 4:
        raise ValueError("expected six patches and four inline-component bridges")
    intervals = []
    minimum_minors = [None, None, None]
    segment_counts = []
    for patch_index, cover in enumerate(covers):
        segments = cover.get("segments", [])
        segment_counts.append(len(segments))
        if not segments:
            raise ValueError("empty metric patch")
        previous_lower = None
        for expected_index, segment in enumerate(segments):
            if segment.get("index") != expected_index:
                raise ValueError("nonconsecutive segment index")
            lower, upper = _segment_interval(segment)
            if not lower < upper:
                raise ValueError("degenerate metric segment")
            # Generator records run from the upper path endpoint downward.
            if previous_lower is not None and upper != previous_lower:
                raise ValueError("gap or overlap inside a metric patch")
            previous_lower = lower
            minors = segment.get("leading_principal_minor_lower_bounds", [])
            if [entry.get("dimension") for entry in minors] != [1, 2, 3]:
                raise ValueError("missing leading principal minor")
            for index, entry in enumerate(minors):
                lower_bound = fraction(entry["real_lower"])
                if lower_bound <= 0 or not entry.get("strictly_positive"):
                    raise ValueError("nonpositive metric minor")
                if minimum_minors[index] is None or lower_bound < minimum_minors[index]:
                    minimum_minors[index] = lower_bound
            if not segment.get("corrected_metric_positive"):
                raise ValueError("metric positivity flag contradicts minors")
            intervals.append((lower, upper, patch_index, expected_index))
    if sum(segment_counts) != 1003:
        raise ValueError("continuation chain does not contain exactly 1003 segments")

    ordered = sorted(intervals)
    merged_lower, merged_upper = ordered[0][0], ordered[0][1]
    if merged_lower != EXPECTED_INTERVAL[0]:
        raise ValueError("wrong continuation lower endpoint")
    for lower, upper, _, _ in ordered[1:]:
        if lower > merged_upper:
            raise ValueError("gap between metric patches")
        merged_upper = max(merged_upper, upper)
    if merged_upper != EXPECTED_INTERVAL[1]:
        raise ValueError("wrong continuation upper endpoint")

    bridge_records = []
    for bridge_index, bridge in enumerate(bridges):
        records = bridge.get("krawczyk", {}).get("component_records", [])
        if len(records) != 83:
            raise ValueError("bridge lacks 83 serialized component inclusions")
        for expected_index, record in enumerate(records):
            if record.get("index") != expected_index:
                raise ValueError("nonconsecutive Krawczyk component")
            image = fraction(record["image_abs_upper_rational"])
            radius = fraction(record["target_radius_rational"])
            if image < 0 or not image < radius or not record.get("strictly_inside"):
                raise ValueError("Krawczyk component is not strictly included")
        summary = bridge["krawczyk"]
        if summary.get("component_count") != 83 or summary.get("strictly_included_component_count") != 83:
            raise ValueError("bridge summary disagrees with component records")
        bridge_records.append(records)
    records = bridge3.get("component_records", [])
    if len(records) != 83:
        raise ValueError("fifth bridge lacks component certificate")
    for expected_index, record in enumerate(records):
        if record.get("index") != expected_index:
            raise ValueError("nonconsecutive fifth-bridge component")
        if not fraction(record["image_abs_upper_rational"]) < fraction(record["target_radius_rational"]):
            raise ValueError("fifth-bridge Krawczyk inclusion failed")
        if not record.get("strictly_inside"):
            raise ValueError("fifth-bridge boolean contradicts rational bounds")

    result = {
        "patch_count": 6,
        "bridge_count": 5,
        "metric_segment_count": sum(segment_counts),
        "segment_counts_by_patch": segment_counts,
        "krawczyk_component_inclusion_count": 5 * 83,
        "certified_depth_interval": [pair(merged_lower), pair(merged_upper)],
        "minimum_leading_principal_minor_lower_bounds": [pair(value) for value in minimum_minors],
        "no_metric_segment_gaps": True,
        "all_metric_minors_strictly_positive": True,
        "all_bridge_components_strictly_included": True,
    }
    if robust is not None:
        if base_root_sha256 is None or material_box_sha256 is None:
            raise ValueError("live anchor hashes are required for the robust chain")
        result["robust_root_chain"] = verify_robust_chain(
            robust, base_root_sha256, material_box_sha256
        )
    return result


def verify_live() -> dict:
    cover_paths = [ROOT / "results" / name for name in COVERS]
    bridge_paths = [ROOT / "results" / name for name in BRIDGES]
    bridge3_path = ROOT / "results" / BRIDGE_3_COMPONENTS
    robust_path = ROOT / "results" / ROBUST_CHAIN
    base_root_path = ROOT / "results" / BASE_ROOT
    material_box_path = ROOT / "results" / MATERIAL_BOX
    nested_box_path = ROOT / "results" / NESTED_MATERIAL_BOX
    provenance_audit_path = ROOT / "results" / TRANSITIVE_PROVENANCE_AUDIT
    stratified_path = ROOT / "results" / STRATIFIED_METRIC_FORMULA
    local_paths = {
        name: ROOT / "results" / name for name in (
            INDEPENDENT_MATERIAL_ROOT, INDEPENDENT_BASELINE,
            INDEPENDENT_CORRECTION, INDEPENDENT_MATERIALITY,
        )
    }
    documents = [json.loads(path.read_text(encoding="utf-8")) for path in cover_paths]
    bridges = [json.loads(path.read_text(encoding="utf-8")) for path in bridge_paths]
    bridge3 = json.loads(bridge3_path.read_text(encoding="utf-8"))
    robust = json.loads(robust_path.read_text(encoding="utf-8"))
    theorem = verify_documents(
        documents, bridges, bridge3, robust,
        sha256(base_root_path), sha256(material_box_path),
    )
    local_hashes = {name: sha256(path) for name, path in local_paths.items()}
    local_hashes[NESTED_MATERIAL_BOX] = sha256(nested_box_path)
    theorem["independent_local_material_box"] = verify_local_material_chain(
        *(json.loads(local_paths[name].read_text(encoding="utf-8")) for name in (
            INDEPENDENT_MATERIAL_ROOT, INDEPENDENT_BASELINE,
            INDEPENDENT_CORRECTION, INDEPENDENT_MATERIALITY,
        )),
        local_hashes,
        json.loads(base_root_path.read_text(encoding="utf-8")),
        json.loads(material_box_path.read_text(encoding="utf-8")),
        json.loads(nested_box_path.read_text(encoding="utf-8")),
    )
    stratified = json.loads(stratified_path.read_text(encoding="utf-8"))
    stratified_status = stratified.get("status", {})
    records = stratified.get("records", [])
    if len(records) != 6 or not all(
        len(record.get("independent_leading_principal_minors", [])) == 3
        and all(fraction(minor["directed_rational_bounds"]["lower"]) > 0
                for minor in record["independent_leading_principal_minors"])
        for record in records
    ):
        raise ValueError("stratified independent metric-formula evidence is incomplete")
    if not (
        stratified_status.get("stratified_metric_formula_replication_passed")
        and stratified_status.get("all_18_primary_minor_lowers_inside_independent_enclosures")
        and not stratified_status.get("all_1003_metric_segments_formula_independently_reimplemented")
        and not stratified_status.get("root_box_existence_independently_reproved")
    ):
        raise ValueError("stratified independent metric-formula status is invalid")
    theorem["stratified_independent_metric_formula_replication"] = {
        "patch_count": len(records),
        "minor_count": sum(len(record["independent_leading_principal_minors"])
                           for record in records),
        "selection_rule": stratified["protocol"]["selection_rule"],
        "all_independent_minors_positive": True,
        "all_primary_minor_lowers_inside_independent_enclosures": True,
        "root_boxes_common_not_independently_reproved": True,
    }
    provenance_audit = json.loads(provenance_audit_path.read_text(encoding="utf-8"))
    provenance_status = provenance_audit.get("status", {})
    publication_package, bundle_ready = publication_bundle_structurally_ready()
    if provenance_status.get("continuation_numerical_inequalities_refuted_by_this_audit"):
        raise ValueError("provenance audit improperly claims numerical refutation")
    theorem["transitive_provenance_audit"] = {
        "visited_node_count": provenance_audit["summary"]["visited_node_count"],
        "resolved_edge_count": provenance_audit["summary"]["resolved_edge_count"],
        "unresolved_edge_count": provenance_audit["summary"]["unresolved_edge_count"],
        "unresolved_proof_edge_count": provenance_audit["summary"]["unresolved_proof_edge_count"],
        "unresolved_non_proof_diagnostic_edge_count": provenance_audit["summary"][
            "unresolved_non_proof_diagnostic_edge_count"
        ],
        "publication_excluded_diagnostic_reference_count": provenance_audit["summary"][
            "publication_excluded_diagnostic_reference_count"
        ],
        "historical_manifest_internal_sha256_record_count": provenance_audit["summary"][
            "historical_manifest_internal_sha256_record_count"
        ],
        "ambiguous_edge_count": provenance_audit["summary"]["ambiguous_edge_count"],
        "all_transitive_sha256_edges_resolve_uniquely": provenance_status[
            "all_transitive_sha256_edges_resolve_uniquely"
        ],
        "all_transitive_proof_sha256_edges_resolve_uniquely": provenance_status[
            "all_transitive_proof_sha256_edges_resolve_uniquely"
        ],
        "publication_diagnostic_exclusion_set_matches_historical_manifest_exactly": provenance_status[
            "publication_diagnostic_exclusion_set_matches_historical_manifest_exactly"
        ],
    }
    paths = (cover_paths + bridge_paths
             + [bridge3_path, robust_path, base_root_path, material_box_path,
                nested_box_path, provenance_audit_path, stratified_path]
             + list(local_paths.values()))
    return {
        "schema_version": 1,
        "theorem": "one base-anchored branch: 7 robust root patches + 7 robust bridges (<0.85 target); 1003 positive-metric segments -> depth in [1.49985, 1.5]; at depth 1.5 an independent 83-real local root and entrywise Fraction verifier prove corrected-metric positivity and >1% correction materiality",
        "verifier_independence": {
            "project_generator_imports": False,
            "third_party_imports": False,
            "arithmetic": "fractions.Fraction over serialized rational outward bounds",
        },
        "proof_data": theorem,
        "provenance": {
            path.relative_to(ROOT).as_posix(): sha256(path) for path in paths
        },
        "publication_bundle": {
            "package_name": PUBLICATION_PACKAGE_NAME,
            "structurally_complete": bundle_ready,
            "license_present": (publication_package / "LICENSE").is_file(),
            "citation_cff_present": (publication_package / "CITATION.cff").is_file(),
            "one_command_verifier_present": (publication_package / "verify_release.py").is_file(),
            "release_manifest_present": (publication_package / "MANIFEST.sha256.json").is_file(),
            "externally_executed_ci_history_present": False,
        },
        "status": {
            "all_1003_serialized_metric_minor_records_exactly_replayed": True,
            "five_83_component_krawczyk_bridges_independently_checked": True,
            "seven_robust_root_patches_independently_checked": True,
            "seven_robust_83_component_bridges_independently_checked": True,
            "robust_root_chain_contraction_target_strictly_met": True,
            "connected_interval_exactly_recomputed": True,
            "publication_continuation_theorem_admitted": True,
            "independent_local_material_root_checked": True,
            "independent_local_corrected_metric_fraction_minors_checked": True,
            "independent_local_correction_materiality_checked": True,
            "publication_local_material_box_theorem_admitted": True,
            "local_material_root_identified_with_continuation_endpoint": True,
            "publication_transitive_provenance_release_admitted": bool(
                provenance_status["publication_transitive_provenance_release_admitted"]
            ),
            "internal_transitive_proof_provenance_complete": bool(
                provenance_status["internal_transitive_proof_dag_complete"]
            ),
            "internal_continuation_proof_bundle_admitted": True,
            "continuation_metric_formula_independently_reimplemented_on_all_1003_segments": False,
            "continuation_metric_formula_independently_reimplemented_on_six_weakest_stratified_segments": True,
            "publication_release_bundle_admitted": bool(
                bundle_ready and provenance_status["publication_transitive_provenance_release_admitted"]
            ),
            "publication_release_readiness_admitted": bool(
                bundle_ready and provenance_status["publication_transitive_provenance_release_admitted"]
            ),
            "global_metric_positivity_certified": False,
        },
    }


def main() -> None:
    payload = verify_live()
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"proof_data": payload["proof_data"], "status": payload["status"]}, indent=2))


if __name__ == "__main__":
    main()
