"""Stdlib-only verifier for the robust local continuation-chain addendum.

It checks provenance, exact rational coverage, and every serialized
componentwise Krawczyk inequality without importing any project generator.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results/2026-08-13_string_5_81_glsm_robust_krawczyk_chain.json"
GENERATOR = ROOT / "models/string_5_81_glsm_robust_krawczyk_patch.py"
BASE_ROOT = ROOT / "results/2026-08-13_string_5_81_glsm_affine_path_krawczyk_ladder.json"
MATERIAL_BOX = ROOT / "results/2026-08-13_string_5_81_glsm_depth15_material_interval_box.json"
EXPECTED_LOWER = Fraction(-3, 20000)
EXPECTED_UPPER = Fraction(0)
EXPECTED_TARGET = Fraction(17, 20)
DIMENSION = 83


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def q(value) -> Fraction:
    if not isinstance(value, dict) or set(value) != {"numerator", "denominator"}:
        raise ValueError("noncanonical rational pair")
    denominator = value["denominator"]
    if not isinstance(value["numerator"], int) or not isinstance(denominator, int):
        raise ValueError("rational pair must contain integers")
    if denominator <= 0:
        raise ValueError("rational denominator must be positive")
    result = Fraction(value["numerator"], denominator)
    if result.numerator != value["numerator"] or result.denominator != denominator:
        raise ValueError("rational pair is not reduced")
    return result


def verify_payload(payload: dict, check_live_provenance: bool = True) -> dict:
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported schema version")
    if check_live_provenance:
        if payload["inputs"]["generator_sha256"] != sha256(GENERATOR):
            raise ValueError("generator provenance mismatch")
        if payload["inputs"]["base_root_sha256"] != sha256(BASE_ROOT):
            raise ValueError("base-root provenance mismatch")
        if payload["inputs"]["material_box_sha256"] != sha256(MATERIAL_BOX):
            raise ValueError("material-box provenance mismatch")
    protocol = payload["protocol"]
    if q(protocol["target_maximum_image_radius_ratio"]) != EXPECTED_TARGET:
        raise ValueError("unexpected contraction target")
    patch_radius = q(protocol["patch_radius"])
    bridge_radius = q(protocol["bridge_radius"])
    if patch_radius <= 0 or bridge_radius <= 0:
        raise ValueError("nonpositive protocol radius")
    patches = payload["patches"]
    bridges = payload["bridges"]
    if len(patches) != 7 or len(bridges) != 7:
        raise ValueError("the robust chain must contain seven patches and seven bridges")
    intervals = []
    maximum_patch_ratio = Fraction(0)
    for patch in patches:
        lo, hi = map(q, patch["protocol"]["t_interval"])
        expected_radius = q(protocol["first_patch_radius"]) if not intervals else patch_radius
        if hi <= lo or (hi - lo) / 2 != expected_radius:
            raise ValueError("invalid patch interval")
        intervals.append((lo, hi))
        radii = list(map(q, patch["krawczyk"]["residual_component_radius_rational_pairs"]))
        uppers = list(map(q, patch["krawczyk"]["component_image_abs_upper_rational_pairs"]))
        if len(radii) != DIMENSION or len(uppers) != DIMENSION:
            raise ValueError("patch component count mismatch")
        if any(radius <= 0 or upper < 0 or upper >= radius
               for upper, radius in zip(uppers, radii)):
            raise ValueError("patch componentwise strict inclusion failed")
        ratio = max(upper / radius for upper, radius in zip(uppers, radii))
        if ratio >= EXPECTED_TARGET:
            raise ValueError("patch contraction target failed")
        maximum_patch_ratio = max(maximum_patch_ratio, ratio)
    intervals.sort()
    if intervals[0][0] != EXPECTED_LOWER or intervals[-1][1] != EXPECTED_UPPER:
        raise ValueError("robust chain endpoints changed")
    if any(right[0] > left[1] for left, right in zip(intervals, intervals[1:])):
        raise ValueError("gap in robust patch cover")

    maximum_bridge_ratio = Fraction(0)
    for bridge in bridges:
        radii = list(map(q, bridge["krawczyk"]["component_target_radius_rational_pairs"]))
        uppers = list(map(q, bridge["krawczyk"]["component_image_abs_upper_rational_pairs"]))
        if len(radii) != DIMENSION or len(uppers) != DIMENSION:
            raise ValueError("bridge component count mismatch")
        if any(radius <= 0 or upper < 0 or upper >= radius
               for upper, radius in zip(uppers, radii)):
            raise ValueError("bridge componentwise strict inclusion failed")
        maximum_bridge_ratio = max(
            maximum_bridge_ratio,
            max(upper / radius for upper, radius in zip(uppers, radii)),
        )
        if not bridge["status"]["same_root_branch_across_local_predictor_change_certified"]:
            raise ValueError("bridge status is false")
    summary_interval = list(map(q, payload["summary"]["covered_t_interval"]))
    if summary_interval != [EXPECTED_LOWER, EXPECTED_UPPER]:
        raise ValueError("summary interval mismatch")
    if not all(payload["status"].values()):
        raise ValueError("artifact status is not fully certified")
    return {
        "patch_count": len(patches),
        "bridge_count": len(bridges),
        "componentwise_patch_inequality_count": len(patches) * DIMENSION,
        "componentwise_bridge_inequality_count": len(bridges) * DIMENSION,
        "covered_t_interval": [str(EXPECTED_LOWER), str(EXPECTED_UPPER)],
        "maximum_patch_ratio_exact": str(maximum_patch_ratio),
        "maximum_patch_ratio": float(maximum_patch_ratio),
        "maximum_bridge_ratio_exact": str(maximum_bridge_ratio),
        "maximum_bridge_ratio": float(maximum_bridge_ratio),
        "independent_serialized_chain_verifier_passed": True,
    }


def main() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    print(json.dumps(verify_payload(payload), indent=2))


if __name__ == "__main__":
    main()
