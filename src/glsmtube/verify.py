"""Dependency-free exact last-mile verifier for the product-tube certificate.

This module verifies rational enclosures produced by the separately distributed
Arb regeneration code. It never treats stored Boolean verdicts or float
diagnostics as proof inputs.
"""

from __future__ import annotations

from fractions import Fraction
import gzip
import hashlib
from itertools import product
import json
import math
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DATA = PACKAGE_ROOT / "data"
ROOT_CERTIFICATE = DATA / "root_certificate.json.gz"
METRIC_CERTIFICATE = DATA / "metric_certificate.json.gz"
TUBE_GRID = DATA / "tube_grid.json"
ROBUST_CHAIN = DATA / "robust_chain.json"
CORNER_TILING = DATA / "corner_tiling.json"
CERTIFICATE = DATA / "certificate.json"
EXACT_GEOMETRY = DATA / "exact_integer_geometry.json"
RAW_ADAPTER = PACKAGE_ROOT / "regeneration" / "data" / "raw_adapter.json"
ATLAS = PACKAGE_ROOT / "regeneration" / "data" / "atlas.json"

PATH_COUNT = 1084
CHILDREN_PER_PARENT = 32
UNKNOWN_DIMENSION = 83
LONGITUDINAL_RADIUS = Fraction(7, 100_000_000)
TRANSVERSE_RADIUS = Fraction(1, 100_000_000)
CHILD_RADIUS = Fraction(1, 200_000_000)
PATH_LOWER = Fraction(-3, 20_000)
PATH_UPPER = Fraction(0)
VERTICES = tuple(product((-1, 1), repeat=5))
ROBUST_TARGET = Fraction(17, 20)


class CertificateError(ValueError):
    """Raised when a proof record fails closed."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def q(value) -> Fraction:
    if not isinstance(value, dict) or set(value) != {"numerator", "denominator"}:
        raise CertificateError("noncanonical rational pair")
    numerator, denominator = value["numerator"], value["denominator"]
    if not isinstance(numerator, int) or isinstance(numerator, bool):
        raise CertificateError("non-integer numerator")
    if not isinstance(denominator, int) or isinstance(denominator, bool) or denominator <= 0:
        raise CertificateError("invalid denominator")
    result = Fraction(numerator, denominator)
    if (result.numerator, result.denominator) != (numerator, denominator):
        raise CertificateError("unreduced rational pair")
    return result


def directed_interval(value):
    if not isinstance(value, dict) or set(value) != {"lower", "upper"}:
        raise CertificateError("invalid directed interval schema")
    lower, upper = q(value["lower"]), q(value["upper"])
    if lower > upper:
        raise CertificateError("reversed directed interval")
    return lower, upper


def exact_integer(value) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CertificateError("non-numeric integer geometry entry")
    if isinstance(value, int):
        return value
    if not value.is_integer() or abs(value) > 2**53:
        raise CertificateError("integer geometry entry would require rounding")
    return int(value)


def verify_exact_geometry(registry, raw_adapter, atlas):
    source = registry.get("source", {})
    if source.get("raw_adapter_sha256") != sha256(RAW_ADAPTER):
        raise CertificateError("raw adapter provenance mismatch")
    if source.get("atlas_sha256") != sha256(ATLAS):
        raise CertificateError("atlas provenance mismatch")
    raw = raw_adapter["toric_sampler_core"]
    charges = [[exact_integer(value) for value in row] for row in raw["glsm_charges"]]
    exponents = [[exact_integer(value) for value in row] for row in raw["exp_aK"]]
    selected = [exact_integer(value) for value in atlas["fan"]["selected_cone_divisor_indices"]]
    geometry = registry.get("geometry", {})
    if geometry.get("glsm_charges") != charges:
        raise CertificateError("exact charge registry mismatch")
    if geometry.get("hypersurface_exponents") != exponents:
        raise CertificateError("exact exponent registry mismatch")
    if geometry.get("selected_cone_divisor_indices_one_based") != selected:
        raise CertificateError("exact chart-index registry mismatch")
    if geometry.get("hypersurface_coefficients_acb_text") != raw["coeff_aK_text"]:
        raise CertificateError("hypersurface coefficient registry mismatch")
    return len(charges), len(exponents)


def verify_positive_minor_intervals(intervals, sequence):
    if not isinstance(intervals, list) or len(intervals) != 3:
        raise CertificateError(f"metric minor count mutation at {sequence}")
    for item in intervals:
        lower, upper = directed_interval(item)
        if not 0 < lower <= upper:
            raise CertificateError(f"metric positivity failure at {sequence}")
    return directed_interval(intervals[2])[0]


def verify_robust_chain(chain):
    if chain.get("schema_version") != 1:
        raise CertificateError("unsupported robust-chain schema")
    protocol = chain["protocol"]
    if q(protocol["target_maximum_image_radius_ratio"]) != ROBUST_TARGET:
        raise CertificateError("robust-chain contraction target mutation")
    patch_radius = q(protocol["patch_radius"])
    bridge_radius = q(protocol["bridge_radius"])
    if patch_radius <= 0 or bridge_radius <= 0:
        raise CertificateError("nonpositive robust-chain protocol radius")
    patches, bridges = chain["patches"], chain["bridges"]
    if len(patches) != 7 or len(bridges) != 7:
        raise CertificateError("unexpected robust-chain size")
    intervals = []
    maximum_patch_ratio = Fraction(0)
    for index, patch in enumerate(patches):
        lower, upper = map(q, patch["protocol"]["t_interval"])
        expected_radius = q(protocol["first_patch_radius"]) if index == 0 else patch_radius
        if upper <= lower or (upper - lower) / 2 != expected_radius:
            raise CertificateError(f"invalid robust patch interval at {index}")
        intervals.append((lower, upper))
        radii = list(map(q, patch["krawczyk"]["residual_component_radius_rational_pairs"]))
        images = list(map(q, patch["krawczyk"]["component_image_abs_upper_rational_pairs"]))
        if len(radii) != UNKNOWN_DIMENSION or len(images) != UNKNOWN_DIMENSION:
            raise CertificateError(f"robust patch dimension mutation at {index}")
        if any(radius <= 0 or image < 0 or image >= radius
               for image, radius in zip(images, radii)):
            raise CertificateError(f"robust patch inclusion failure at {index}")
        ratio = max(image / radius for image, radius in zip(images, radii))
        if ratio >= ROBUST_TARGET:
            raise CertificateError(f"robust patch contraction failure at {index}")
        maximum_patch_ratio = max(maximum_patch_ratio, ratio)
    intervals.sort()
    if intervals[0][0] != PATH_LOWER or intervals[-1][1] != PATH_UPPER:
        raise CertificateError("robust-chain endpoint mutation")
    if any(right[0] > left[1] for left, right in zip(intervals, intervals[1:])):
        raise CertificateError("gap in robust patch cover")

    maximum_bridge_ratio = Fraction(0)
    for index, bridge in enumerate(bridges):
        radii = list(map(q, bridge["krawczyk"]["component_target_radius_rational_pairs"]))
        images = list(map(q, bridge["krawczyk"]["component_image_abs_upper_rational_pairs"]))
        if len(radii) != UNKNOWN_DIMENSION or len(images) != UNKNOWN_DIMENSION:
            raise CertificateError(f"robust bridge dimension mutation at {index}")
        if any(radius <= 0 or image < 0 or image >= radius
               for image, radius in zip(images, radii)):
            raise CertificateError(f"robust bridge inclusion failure at {index}")
        if bridge.get("status", {}).get(
                "same_root_branch_across_local_predictor_change_certified") is not True:
            raise CertificateError(f"robust bridge branch mutation at {index}")
        maximum_bridge_ratio = max(
            maximum_bridge_ratio,
            max(image / radius for image, radius in zip(images, radii)),
        )
    if list(map(q, chain["summary"]["covered_t_interval"])) != [PATH_LOWER, PATH_UPPER]:
        raise CertificateError("robust-chain summary mutation")
    return maximum_patch_ratio, maximum_bridge_ratio


def derive_centres(patch, t):
    local = patch["local_affine_coordinates"]
    slopes = list(map(q, local["slope_rational_pairs"]))
    intercepts = list(map(q, local["intercept_rational_pairs"]))
    if len(slopes) != UNKNOWN_DIMENSION or len(intercepts) != UNKNOWN_DIMENSION:
        raise CertificateError("robust predictor dimension mutation")
    return [intercept + slope * t for intercept, slope in zip(intercepts, slopes)]


def verify_tiling(payload):
    exact = payload["exact_tiling"]
    if exact["dimension"] != 5 or exact["corner_alphabet"] != [-1, 1]:
        raise CertificateError("corner geometry mutation")
    boxes = exact["corner_boxes"]
    if len(boxes) != CHILDREN_PER_PARENT:
        raise CertificateError("corner child count mutation")
    for index, (box, vertex) in enumerate(zip(boxes, VERTICES)):
        if box["sequence"] != index or box["offset_slots"] != list(vertex):
            raise CertificateError(f"corner ordering mutation at {index}")
        centres = list(map(q, box["center_offsets"]))
        widths = list(map(q, box["half_widths"]))
        expected_centres = [slot * CHILD_RADIUS for slot in vertex]
        if centres != expected_centres or widths != [CHILD_RADIUS] * 5:
            raise CertificateError(f"corner interval mutation at {index}")
        for centre, width, slot in zip(centres, widths, vertex):
            expected = (-TRANSVERSE_RADIUS, 0) if slot == -1 else (0, TRANSVERSE_RADIUS)
            if (centre - width, centre + width) != expected:
                raise CertificateError(f"corner box fails to tile at {index}")
    return len(boxes)


def verify_path_grid(tube):
    protocol, records = tube["protocol"], tube["records"]
    if protocol["center_count"] != PATH_COUNT or len(records) != PATH_COUNT:
        raise CertificateError("path count mutation")
    if q(protocol["t_lower"]) != PATH_LOWER or q(protocol["t_upper"]) != PATH_UPPER:
        raise CertificateError("path endpoint mutation")
    spacing = q(protocol["center_spacing"])
    if 2 * LONGITUDINAL_RADIUS < spacing:
        raise CertificateError("longitudinal gap")
    intervals = []
    for index, record in enumerate(records):
        expected = PATH_LOWER + index * spacing
        if q(record["t_center"]) != expected:
            raise CertificateError(f"path ordering mutation at {index}")
        if q(record["depth_center"]) != Fraction(3, 2) + expected:
            raise CertificateError(f"depth mutation at {index}")
        intervals.append((max(PATH_LOWER, expected - LONGITUDINAL_RADIUS),
                          min(PATH_UPPER, expected + LONGITUDINAL_RADIUS)))
    if q(records[-1]["t_center"]) != PATH_UPPER:
        raise CertificateError("path does not reach upper endpoint")
    minimum_overlap = None
    for index, (left, right) in enumerate(zip(intervals, intervals[1:])):
        width = min(left[1], right[1]) - max(left[0], right[0])
        if width <= 0:
            raise CertificateError(f"longitudinal gluing gap after {index}")
        minimum_overlap = width if minimum_overlap is None else min(minimum_overlap, width)
    return intervals, minimum_overlap


def verify_parent(record, index, tube_record, chain):
    if record.get("path_index") != index:
        raise CertificateError(f"parent ordering mutation at {index}")
    if record.get("t_center") != tube_record["t_center"] or record.get("depth_center") != tube_record["depth_center"]:
        raise CertificateError(f"parent coordinate mutation at {index}")
    t = q(tube_record["t_center"])
    patch_index = record.get("robust_patch_index")
    if not isinstance(patch_index, int) or not 0 <= patch_index < len(chain["patches"]):
        raise CertificateError(f"invalid robust patch at {index}")
    patch = chain["patches"][patch_index]
    patch_lower, patch_upper = map(q, patch["protocol"]["t_interval"])
    if not patch_lower <= t <= patch_upper:
        raise CertificateError(f"parent anchor escapes robust patch at {index}")
    centres = list(map(q, record["unknown_center_rational_pairs"]))
    if centres != derive_centres(patch, t):
        raise CertificateError(f"parent predictor mutation at {index}")
    radii = list(map(q, record["final_unknown_radius_rational_pairs"]))
    if len(radii) != UNKNOWN_DIMENSION or any(radius <= 0 for radius in radii):
        raise CertificateError(f"parent radius mutation at {index}")

    derivative = record["hypersurface_derivative_directed_rational_bounds"]
    derivative_real = directed_interval(derivative["real"])
    derivative_imag = directed_interval(derivative["imag"])
    if derivative_real[0] <= 0 <= derivative_real[1] and derivative_imag[0] <= 0 <= derivative_imag[1]:
        raise CertificateError(f"hypersurface derivative may contain zero at {index}")

    newton = record["hypersurface_newton_image_directed_rational_bounds"]
    bounds = [newton[component] for component in ("real", "imag")]
    dterm = record["dterm_krawczyk_image_directed_rational_bounds"]
    if not isinstance(dterm, list) or len(dterm) != 81:
        raise CertificateError(f"D-term image dimension mutation at {index}")
    bounds.extend(dterm)
    margins = []
    for component, (centre, radius, directed) in enumerate(zip(centres, radii, bounds)):
        lower, upper = directed_interval(directed)
        lower_margin = lower - (centre - radius)
        upper_margin = (centre + radius) - upper
        if lower_margin <= 0 or upper_margin <= 0:
            raise CertificateError(f"strict inclusion failure at {index}:{component}")
        margins.extend((lower_margin, upper_margin))

    robust_radii = list(map(q, patch["krawczyk"]["residual_component_radius_rational_pairs"]))
    if len(robust_radii) != UNKNOWN_DIMENSION or any(old >= new for old, new in zip(robust_radii, radii)):
        raise CertificateError(f"robust anchor is not strictly nested at {index}")
    return min(margins), max(old / new for old, new in zip(robust_radii, radii))


def _stream_metric_records(path: Path):
    """Stream the large pretty-printed metric JSON without loading it in memory."""
    marker = '  "records": ['
    handle = gzip.open(path, "rt", encoding="utf-8")
    prefix = ""
    while marker not in prefix:
        block = handle.read(65536)
        if not block:
            handle.close()
            raise CertificateError("metric record marker missing")
        prefix += block
    before, buffer = prefix.split(marker, 1)
    header = json.loads(before + '  "records": []\n}')
    decoder = json.JSONDecoder()
    trailer = {}

    def records():
        nonlocal buffer
        try:
            while True:
                buffer = buffer.lstrip()
                if buffer.startswith("]"):
                    suffix = buffer[1:] + handle.read()
                    trailer["value"] = json.loads('{"records": []' + suffix)
                    return
                if buffer.startswith(","):
                    buffer = buffer[1:]
                    continue
                try:
                    value, used = decoder.raw_decode(buffer)
                except json.JSONDecodeError:
                    block = handle.read(65536)
                    if not block:
                        raise CertificateError("truncated metric record array")
                    buffer += block
                    continue
                buffer = buffer[used:]
                yield value
        finally:
            handle.close()

    return header, records(), trailer


def verify_metric(metric_path, tube):
    header, records, trailer_holder = _stream_metric_records(metric_path)
    protocol = header["protocol"]
    total = PATH_COUNT * CHILDREN_PER_PARENT
    if protocol["path_center_count"] != PATH_COUNT or protocol["record_count"] != total:
        raise CertificateError("metric protocol count mutation")
    if q(protocol["longitudinal_radius"]) != LONGITUDINAL_RADIUS or q(protocol["transverse_radius"]) != CHILD_RADIUS:
        raise CertificateError("metric protocol geometry mutation")
    minimum_minors = [None, None, None]
    maximum_q = 0.0
    count = 0
    for sequence, record in enumerate(records):
        if sequence >= total:
            raise CertificateError("excess metric records")
        parent, child = divmod(sequence, CHILDREN_PER_PARENT)
        if record.get("sequence") != sequence or record.get("path_index") != parent:
            raise CertificateError(f"metric sequence mutation at {sequence}")
        if record.get("offset_slots") != list(VERTICES[child]):
            raise CertificateError(f"metric corner mutation at {sequence}")
        tube_record = tube["records"][parent]
        if record.get("t_center") != tube_record["t_center"] or record.get("depth_center") != tube_record["depth_center"]:
            raise CertificateError(f"metric path mutation at {sequence}")
        intervals = record.get("minor_directed_rational_bounds")
        verify_positive_minor_intervals(intervals, sequence)
        for minor_index, interval in enumerate(intervals):
            lower = directed_interval(interval)[0]
            current = minimum_minors[minor_index]
            minimum_minors[minor_index] = lower if current is None else min(current, lower)
        root_q = record.get("root_q")
        if (not isinstance(root_q, (int, float)) or isinstance(root_q, bool)
                or not math.isfinite(root_q)):
            raise CertificateError(f"invalid diagnostic q at {sequence}")
        maximum_q = max(maximum_q, float(root_q))
        count += 1
    if count != total:
        raise CertificateError(f"metric record count {count} != {total}")
    trailer = trailer_holder.get("value")
    if trailer is None:
        raise CertificateError("metric trailer missing")
    summary = trailer["summary"]
    if summary["root_pass_count"] != total or summary["metric_pass_count"] != total:
        raise CertificateError("metric summary count mutation")
    if Fraction(summary["minimum_third_minor_lower_exact"]) != minimum_minors[2]:
        raise CertificateError("metric minimum mutation")
    return count, minimum_minors, maximum_q


def verify_package():
    from glsmtube.endpoint_verify import verify_payload as verify_endpoint_payload

    certificate = load_json(CERTIFICATE)
    tube = load_json(TUBE_GRID)
    chain = load_json(ROBUST_CHAIN)
    tiling = load_json(CORNER_TILING)
    roots = load_json(ROOT_CERTIFICATE)
    geometry_shape = verify_exact_geometry(
        load_json(EXACT_GEOMETRY), load_json(RAW_ADAPTER), load_json(ATLAS)
    )
    robust_patch_ratio, robust_bridge_ratio = verify_robust_chain(chain)
    endpoint = verify_endpoint_payload(load_json(
        DATA / "endpoint_validation" /
        "2026-08-13_string_5_81_glsm_depth15_open_neighborhood_certificate.json"
    ))

    if roots.get("schema_version") != 1:
        raise CertificateError("unsupported root certificate schema")
    protocol = roots["protocol"]
    if protocol["parent_root_box_count"] != PATH_COUNT or protocol["unknown_dimension"] != UNKNOWN_DIMENSION:
        raise CertificateError("root protocol mutation")
    if q(protocol["longitudinal_radius"]) != LONGITUDINAL_RADIUS or q(protocol["full_transverse_radius"]) != TRANSVERSE_RADIUS:
        raise CertificateError("root domain mutation")

    _, minimum_overlap = verify_path_grid(tube)
    child_templates = verify_tiling(tiling)
    records = roots["records"]
    if len(records) != PATH_COUNT:
        raise CertificateError("parent record count mutation")
    minimum_margin = None
    maximum_anchor_ratio = Fraction(0)
    for index, (record, tube_record) in enumerate(zip(records, tube["records"])):
        margin, anchor_ratio = verify_parent(record, index, tube_record, chain)
        minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)
        maximum_anchor_ratio = max(maximum_anchor_ratio, anchor_ratio)

    metric_count, minimum_minors, maximum_q = verify_metric(METRIC_CERTIFICATE, tube)

    theorem = certificate["theorem"]
    expected = {
        "parent_exact_root_box_count": PATH_COUNT,
        "metric_child_box_count": metric_count,
        "exact_componentwise_root_inclusion_count": PATH_COUNT * UNKNOWN_DIMENSION,
        "longitudinal_gluing_edge_count": PATH_COUNT - 1,
        "minimum_componentwise_root_margin_exact": str(minimum_margin),
        "minimum_metric_first_minor_lower_exact": str(minimum_minors[0]),
        "minimum_metric_second_minor_lower_exact": str(minimum_minors[1]),
        "minimum_metric_third_minor_lower_exact": str(minimum_minors[2]),
        "minimum_longitudinal_overlap_width_exact": str(minimum_overlap),
        "maximum_robust_anchor_to_parent_radius_ratio_exact": str(maximum_anchor_ratio),
    }
    for key, value in expected.items():
        if theorem.get(key) != value:
            raise CertificateError(f"summary mismatch: {key}")
    scope = certificate["scope"]
    if any(scope[name] for name in (
        "global_cy_positivity_certified", "ricci_flatness_certified",
        "physical_vacuum_certified"
    )):
        raise CertificateError("scope overclaim")
    return {
        "verified": True,
        "parent_root_boxes": PATH_COUNT,
        "metric_child_boxes": metric_count,
        "corner_templates_per_parent": child_templates,
        "exact_root_components": PATH_COUNT * UNKNOWN_DIMENSION,
        "longitudinal_gluing_edges": PATH_COUNT - 1,
        "minimum_root_margin_exact": str(minimum_margin),
        "minimum_metric_minors_exact": list(map(str, minimum_minors)),
        "minimum_longitudinal_overlap_exact": str(minimum_overlap),
        "maximum_anchor_ratio_exact": str(maximum_anchor_ratio),
        "maximum_robust_patch_ratio_exact": str(robust_patch_ratio),
        "maximum_robust_bridge_ratio_exact": str(robust_bridge_ratio),
        "maximum_root_q_float_diagnostic_only": maximum_q,
        "exact_integer_geometry_shapes": geometry_shape,
        "independent_endpoint_validation": endpoint,
        "quick_verifier_scope": "exact last-mile replay; Arb formula regeneration is separate",
    }


def main():
    print(json.dumps(verify_package(), indent=2))


if __name__ == "__main__":
    main()
