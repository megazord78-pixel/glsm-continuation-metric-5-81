"""Independent hypersurface + D-term root proof for the weakest tube box."""

from __future__ import annotations

import csv
from fractions import Fraction
import hashlib
import json
import re
from pathlib import Path

from regeneration.models.string_5_81_glsm_stratified_metric_formula_verifier import (
    exact,
    exact_integer,
    fraction,
    pair,
    parse_acb,
    rational_bounds,
)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
GENERATOR = Path(__file__).resolve()
INPUTS = RESULTS / "2026-08-13_string_5_81_glsm_transverse_tube_independent_inputs.json"
ADAPTER = RESULTS / "2026-08-12_string_5_81_cymetric_toric_adapter.json"
ATLAS = RESULTS / "2026-08-12_string_5_81_d3_local_atlas_witness.json"
DECOMPOSITION = RESULTS / "2026-08-12_string_5_81_nef_decomposition_search.json"
KAHLER = ROOT / "research/source_cache/2107.09064/anc/paper_data/5-81-3213/corrected_kahler_param.dat"
RESULT = RESULTS / "2026-08-13_string_5_81_glsm_transverse_tube_independent_root_verifier.json"
PRECISION_DIGITS = 100


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_fi(arb):
    basis = json.loads(DECOMPOSITION.read_text(encoding="utf-8"))["basis"]["source_to_cytools"]
    with KAHLER.open(newline="", encoding="utf-8") as handle:
        text = [value for row in csv.reader(handle) for value in row if value.strip()]
    target = [arb(value) for value in text]
    if len(target) != 81:
        raise ValueError("expected 81 exact Kahler target entries")
    return [
        sum((arb(str(value)) * target_value for value, target_value in zip(row, target)), arb(0))
        for row in basis
    ]


def verify_record(record, adapter, atlas, arb, acb, arb_mat, fmpq):
    selected = [int(value) - 1 for value in atlas["fan"]["selected_cone_divisor_indices"]]
    free_indices = [selected[0], selected[1], selected[3]]
    eliminated_index = selected[2]
    exponents = [[exact_integer(x) for x in row] for row in adapter["exp_aK"]]
    coefficients = [parse_acb(text, arb, acb) for text in adapter["coeff_aK_text"]]
    charges = [[exact_integer(x) for x in row] for row in adapter["glsm_charges"]]

    free_centers = [exact(fraction(x), arb, fmpq) for x in record["free_center_rational_pairs"]]
    free_radii = [exact(fraction(x), arb, fmpq) for x in record["free_radius_rational_pairs"]]
    free = [acb(arb(free_centers[i], free_radii[i]),
                arb(free_centers[i + 3], free_radii[i + 3])) for i in range(3)]
    unknown_centers = [exact(fraction(x), arb, fmpq) for x in record["unknown_center_rational_pairs"]]
    uniqueness_radii = [exact(fraction(x), arb, fmpq) for x in record["unknown_radius_rational_pairs"]]
    eliminated_center = acb(unknown_centers[0], unknown_centers[1])
    eliminated_box = acb(
        arb(unknown_centers[0], uniqueness_radii[0]),
        arb(unknown_centers[1], uniqueness_radii[1]),
    )

    def polynomial(eliminated, derivative=False):
        total = acb(0)
        for coefficient, exponent in zip(coefficients, exponents):
            argument = sum((exponent[index] * free[local]
                            for local, index in enumerate(free_indices)), acb(0))
            argument += exponent[eliminated_index] * eliminated
            term = coefficient * argument.exp()
            total += term * (exponent[eliminated_index] if derivative else 1)
        return total

    derivative = polynomial(eliminated_box, derivative=True)
    if derivative.contains(0):
        raise ValueError("independent hypersurface derivative contains zero")
    newton_image = eliminated_center - polynomial(eliminated_center) / derivative
    hypersurface_included = eliminated_box.contains_interior(newton_image)
    newton_bounds = {
        "real": rational_bounds(newton_image.real, arb, fmpq),
        "imag": rational_bounds(newton_image.imag, arb, fmpq),
    }

    dterm_centers = unknown_centers[2:]
    dterm_radii = uniqueness_radii[2:]
    dterm_box = [arb(center, radius) for center, radius in zip(dterm_centers, dterm_radii)]
    fi = exact_fi(arb)

    log_abs_squared = [arb(0) for _ in range(85)]
    for local, index in enumerate(free_indices):
        log_abs_squared[index] = 2 * free[local].real
    log_abs_squared[eliminated_index] = 2 * newton_image.real

    def weights(values):
        return [
            (log_abs_squared[column] + sum((
                arb(charges[row][column]) * values[row] for row in range(81)
            ), arb(0))).exp()
            for column in range(85)
        ]

    def jacobian(weight_values):
        return arb_mat([[
            sum((arb(charges[row][column] * charges[col][column]) * weight_values[column]
                 for column in range(85)), arb(0))
            for col in range(81)] for row in range(81)
        ])

    center_weights = weights(dterm_centers)
    midpoint_jacobian = jacobian([value.mid() for value in center_weights])
    inverse = midpoint_jacobian.inv()
    preconditioner = arb_mat([[inverse[i, j].mid() for j in range(81)] for i in range(81)])
    residual = arb_mat([[
        sum((arb(charges[row][column]) * center_weights[column]
             for column in range(85)), arb(0)) - fi[row]
    ] for row in range(81)])
    interval_jacobian = jacobian(weights(dterm_box))
    identity = arb_mat(81, 81)
    for i in range(81):
        identity[i, i] = 1
    center_vector = arb_mat([[value] for value in dterm_centers])
    displacement = arb_mat([[arb(0, radius)] for radius in dterm_radii])
    image = center_vector - preconditioner * residual + (identity - preconditioner * interval_jacobian) * displacement
    inclusions = [dterm_box[i].contains_interior(image[i, 0]) for i in range(81)]
    ratios = [float(abs(image[i, 0] - dterm_centers[i]).upper()) / float(dterm_radii[i]) for i in range(81)]
    image_bounds = [rational_bounds(image[i, 0], arb, fmpq) for i in range(81)]
    return {
        "tube_record_index": record["tube_record_index"],
        "hypersurface": {
            "derivative_excludes_zero": not derivative.contains(0),
            "newton_image_strictly_included": hypersurface_included,
            "newton_image_directed_rational_bounds": newton_bounds,
        },
        "dterm": {
            "component_count": 81,
            "strictly_included_component_count": sum(inclusions),
            "maximum_image_radius_ratio": max(ratios),
            "image_directed_rational_bounds": image_bounds,
        },
        "root_existence_independently_reproved": hypersurface_included and all(inclusions),
    }


def refined_verify_record(record, adapter, atlas, arb, acb, arb_mat, fmpq,
                          refinement_steps=3, safety_factor=Fraction(6, 5)):
    """Independently refine D-term radii from a uniform uniqueness box."""
    working = json.loads(json.dumps(record))
    history = []
    for step in range(refinement_steps + 1):
        result = verify_record(working, adapter, atlas, arb, acb, arb_mat, fmpq)
        history.append({
            "step": step,
            "included_component_count": result["dterm"]["strictly_included_component_count"],
            "maximum_image_radius_ratio": result["dterm"]["maximum_image_radius_ratio"],
        })
        if not result["root_existence_independently_reproved"] or step == refinement_steps:
            result["independent_refinement_history"] = history
            result["final_unknown_radius_rational_pairs"] = working["unknown_radius_rational_pairs"]
            return result
        centers = list(map(fraction, working["unknown_center_rational_pairs"]))
        radii = list(map(fraction, working["unknown_radius_rational_pairs"]))
        next_radii = radii[:2]
        for center, bounds in zip(centers[2:], result["dterm"]["image_directed_rational_bounds"]):
            lower, upper = fraction(bounds["lower"]), fraction(bounds["upper"])
            radius = safety_factor * max(abs(lower - center), abs(upper - center))
            next_radii.append(radius)
        working["unknown_radius_rational_pairs"] = [pair(value) for value in next_radii]
    raise AssertionError("unreachable")


def compute():
    from flint import acb, arb, arb_mat, ctx, fmpq
    ctx.dps = PRECISION_DIGITS
    frozen = json.loads(INPUTS.read_text(encoding="utf-8"))
    weakest = frozen["protocol"]["global_weakest_record_index"]
    adapter = json.loads(ADAPTER.read_text(encoding="utf-8"))["toric_sampler_core"]
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    records = [verify_record(record, adapter, atlas, arb, acb, arb_mat, fmpq)
               for record in frozen["records"]]
    passed = all(item["root_existence_independently_reproved"] for item in records)
    return {
        "schema_version": 1,
        "inputs": {
            "generator_sha256": sha256(GENERATOR),
            "frozen_tube_inputs_sha256": sha256(INPUTS),
            "adapter_sha256": sha256(ADAPTER),
            "atlas_sha256": sha256(ATLAS),
            "decomposition_sha256": sha256(DECOMPOSITION),
            "corrected_kahler_parameter_sha256": sha256(KAHLER),
        },
        "protocol": {
            "tube_record_indices": [item["tube_record_index"] for item in records],
            "global_weakest_record_index": weakest,
            "imports_primary_hypersurface_or_moment_map_generator": False,
            "sequential_hypersurface_then_dterm_krawczyk": True,
            "raw_coefficients_charges_fi_rebuilt": True,
            "arb_precision_decimal_digits": PRECISION_DIGITS,
        },
        "records": records,
        "status": {
            "weakest_tube_box_root_existence_independently_reproved": next(
                item["root_existence_independently_reproved"] for item in records
                if item["tube_record_index"] == weakest
            ),
            "all_eight_stratified_root_boxes_independently_reproved": len(records) == 8 and passed,
            "all_1084_root_boxes_independently_reproved": False,
            "global_corrected_metric_positivity_certified": False,
        },
    }


def main():
    payload = compute()
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    compact = [{
        "tube_record_index": item["tube_record_index"],
        "hypersurface_included": item["hypersurface"]["newton_image_strictly_included"],
        "dterm_included": item["dterm"]["strictly_included_component_count"],
        "dterm_q": item["dterm"]["maximum_image_radius_ratio"],
        "passed": item["root_existence_independently_reproved"],
    } for item in payload["records"]]
    print(json.dumps({"protocol": payload["protocol"], "records": compact, "status": payload["status"]}, indent=2))


if __name__ == "__main__":
    main()
