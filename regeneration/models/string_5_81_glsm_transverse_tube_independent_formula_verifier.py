"""Independent interval metric reconstruction on frozen six-real tube boxes.

This module does not import the primary hypersurface, moment-map, baseline-
metric, corrected-metric, or tube generators.  It reuses only generic Jet and
verified-inverse helpers from the older independent formula implementation.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

from regeneration.models.string_5_81_glsm_stratified_metric_formula_verifier import (
    FREE_DIMENSION,
    UNKNOWN_DIMENSION,
    Jet,
    exact,
    exact_integer,
    fraction,
    pair,
    parse_acb,
    rational_bounds,
    verified_inverse,
)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
GENERATOR = Path(__file__).resolve()
INPUTS = RESULTS / "2026-08-13_string_5_81_glsm_transverse_tube_independent_inputs.json"
ADAPTER = RESULTS / "2026-08-12_string_5_81_cymetric_toric_adapter.json"
ATLAS = RESULTS / "2026-08-12_string_5_81_d3_local_atlas_witness.json"
PARAMETERS = RESULTS / "2026-08-13_string_5_81_glsm_boundary_safe_exact_parameters.json"
RESULT = RESULTS / "2026-08-13_string_5_81_glsm_transverse_tube_independent_formula_verifier.json"
PRECISION_DIGITS = 100


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def raw_inputs(arb, acb):
    raw = json.loads(ADAPTER.read_text(encoding="utf-8"))["toric_sampler_core"]
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    selected = [int(x) - 1 for x in atlas["fan"]["selected_cone_divisor_indices"]]
    return {
        "charges": [[exact_integer(x) for x in row] for row in raw["glsm_charges"]],
        "exponents": [[exact_integer(x) for x in row] for row in raw["exp_aK"]],
        "coefficients": [parse_acb(x, arb, acb) for x in raw["coeff_aK_text"]],
        "free_indices": [selected[0], selected[1], selected[3]],
        "eliminated": selected[2],
    }


def finite_abs_upper(value):
    bound = abs(value).upper()
    return float(bound) if bound.is_finite() else None


def maximum_finite_or_infinity(values):
    bounds = [finite_abs_upper(value) for value in values]
    return float("inf") if any(value is None for value in bounds) else max(bounds)


def exact_interval(bounds, arb, fmpq):
    lower = fraction(bounds["lower"])
    upper = fraction(bounds["upper"])
    if not lower <= upper:
        raise ValueError("reversed directed interval")
    midpoint = (lower + upper) / 2
    radius = (upper - lower) / 2
    return arb(exact(midpoint, arb, fmpq), exact(radius, arb, fmpq)), exact(radius, arb, fmpq)


def verified_solve(matrix, rhs, component_radii, arb, arb_mat, order=30):
    """Rigorous scaled Neumann solve without constructing a full inverse."""
    n = matrix.nrows()
    midpoint = arb_mat([[matrix[i, j].mid() for j in range(n)] for i in range(n)])
    preconditioner = midpoint.inv()
    identity = arb_mat(n, n)
    for i in range(n):
        identity[i, i] = 1
    remainder = identity - preconditioner * matrix
    scaled = arb_mat([[
        remainder[i, j] * component_radii[j] / component_radii[i]
        for j in range(n)] for i in range(n)
    ])
    rho = max(sum((abs(scaled[i, j]) for j in range(n)), arb(0)) for i in range(n))
    if not rho < 1:
        raise ValueError(f"independent scaled solve norm is not below one: {rho}")
    initial_unscaled = preconditioner * rhs
    initial = arb_mat([[
        initial_unscaled[i, column] / component_radii[i]
        for column in range(rhs.ncols())
    ] for i in range(n)])
    initial_norm = max(
        sum((abs(initial[i, column]) for column in range(rhs.ncols())), arb(0))
        for i in range(n)
    )
    total = initial
    term = initial
    for _ in range(order):
        term = scaled * term
        total = total + term
    rho_upper = arb(rho.upper())
    tail = rho_upper ** (order + 1) / (1 - rho_upper) * initial_norm
    solved = arb_mat([[
        component_radii[i] * (total[i, column] + arb(0, tail))
        for column in range(rhs.ncols())
    ] for i in range(n)])
    return solved, rho, initial_norm, tail


def verify_box(record, arb, acb, acb_mat, arb_mat, fmpq):
    raw = raw_inputs(arb, acb)
    free_center = [exact(fraction(x), arb, fmpq) for x in record["free_center_rational_pairs"]]
    free_radius = [exact(fraction(x), arb, fmpq) for x in record["free_radius_rational_pairs"]]
    free = [acb(arb(free_center[i], free_radius[i]),
                arb(free_center[i + 3], free_radius[i + 3])) for i in range(3)]
    eliminated_bounds = record["primary_hypersurface_newton_image_directed_rational_bounds"]
    unknown_with_radii = [
        exact_interval(eliminated_bounds[component], arb, fmpq)
        for component in ("real", "imag")
    ]
    unknown_with_radii.extend(
        exact_interval(bounds, arb, fmpq)
        for bounds in record["primary_dterm_krawczyk_image_directed_rational_bounds"]
    )
    unknown = [item[0] for item in unknown_with_radii]
    radii = [item[1] for item in unknown_with_radii]
    if len(unknown) != UNKNOWN_DIMENSION or not all(radius > 0 for radius in radii):
        raise ValueError("invalid certified Krawczyk-image box")

    logs = [acb(0) for _ in range(85)]
    for local, index in enumerate(raw["free_indices"]):
        logs[index] = free[local]
    logs[raw["eliminated"]] = acb(unknown[0], unknown[1])
    terms = []
    for coefficient, exponent in zip(raw["coefficients"], raw["exponents"]):
        argument = sum((exponent[i] * logs[i] for i in range(85)), acb(0))
        terms.append(coefficient * argument.exp())
    weights = []
    for column in range(85):
        log_abs = arb(0)
        for local, index in enumerate(raw["free_indices"]):
            if column == index:
                log_abs += 2 * free[local].real
        if column == raw["eliminated"]:
            log_abs += 2 * unknown[0]
        log_abs += sum((raw["charges"][row][column] * unknown[row + 2]
                        for row in range(81)), arb(0))
        weights.append(log_abs.exp())

    poly_x = [[acb(0) for _ in range(FREE_DIMENSION)] for _ in terms]
    poly_u = [[acb(0) for _ in range(UNKNOWN_DIMENSION)] for _ in terms]
    weight_x = [[arb(0) for _ in range(FREE_DIMENSION)] for _ in range(85)]
    weight_u = [[arb(0) for _ in range(UNKNOWN_DIMENSION)] for _ in range(85)]
    for m, exponent in enumerate(raw["exponents"]):
        for local, index in enumerate(raw["free_indices"]):
            poly_x[m][local] = acb(exponent[index])
            poly_x[m][local + 3] = acb(0, exponent[index])
        poly_u[m][0] = acb(exponent[raw["eliminated"]])
        poly_u[m][1] = acb(0, exponent[raw["eliminated"]])
    for column in range(85):
        for local, index in enumerate(raw["free_indices"]):
            if column == index:
                weight_x[column][local] = arb(2)
        if column == raw["eliminated"]:
            weight_u[column][0] = arb(2)
        for row in range(81):
            weight_u[column][row + 2] = arb(raw["charges"][row][column])

    jacobian = arb_mat(UNKNOWN_DIMENSION, UNKNOWN_DIMENSION)
    free_derivative = arb_mat(UNKNOWN_DIMENSION, FREE_DIMENSION)
    for variable in range(UNKNOWN_DIMENSION):
        derivative = sum((terms[m] * poly_u[m][variable] for m in range(len(terms))), acb(0))
        jacobian[0, variable], jacobian[1, variable] = derivative.real, derivative.imag
        for row in range(81):
            jacobian[row + 2, variable] = sum((
                raw["charges"][row][column] * weights[column] * weight_u[column][variable]
                for column in range(85)), arb(0))
    for variable in range(FREE_DIMENSION):
        derivative = sum((terms[m] * poly_x[m][variable] for m in range(len(terms))), acb(0))
        free_derivative[0, variable], free_derivative[1, variable] = derivative.real, derivative.imag
        for row in range(81):
            free_derivative[row + 2, variable] = sum((
                raw["charges"][row][column] * weights[column] * weight_x[column][variable]
                for column in range(85)), arb(0))
    first_unknown, joint_rho, joint_initial_norm, joint_tail = verified_solve(
        jacobian, -free_derivative, radii, arb, arb_mat
    )
    nonfinite_first = [
        (i, j, str(first_unknown[i, j]))
        for i in range(UNKNOWN_DIMENSION) for j in range(FREE_DIMENSION)
        if not first_unknown[i, j].is_finite()
    ]
    first_unknown_max = maximum_finite_or_infinity(
        first_unknown[i, j] for i in range(UNKNOWN_DIMENSION) for j in range(FREE_DIMENSION)
    )
    poly_direction = [[poly_x[m][i] + sum((poly_u[m][u] * first_unknown[u, i]
                       for u in range(UNKNOWN_DIMENSION)), acb(0))
                       for i in range(FREE_DIMENSION)] for m in range(len(terms))]
    weight_direction = [[weight_x[c][i] + sum((weight_u[c][u] * first_unknown[u, i]
                         for u in range(UNKNOWN_DIMENSION)), arb(0))
                         for i in range(FREE_DIMENSION)] for c in range(85)]
    second = [[None for _ in range(FREE_DIMENSION)] for _ in range(FREE_DIMENSION)]
    for i in range(FREE_DIMENSION):
        for j in range(i, FREE_DIMENSION):
            rhs = arb_mat(UNKNOWN_DIMENSION, 1)
            value = sum((terms[m] * poly_direction[m][i] * poly_direction[m][j]
                         for m in range(len(terms))), acb(0))
            rhs[0, 0], rhs[1, 0] = value.real, value.imag
            for row in range(81):
                rhs[row + 2, 0] = sum((raw["charges"][row][c] * weights[c]
                    * weight_direction[c][i] * weight_direction[c][j] for c in range(85)), arb(0))
            solved, _, _, _ = verified_solve(jacobian, -rhs, radii, arb, arb_mat)
            second[i][j] = second[j][i] = solved
    second_unknown_max = maximum_finite_or_infinity(
        second[i][j][u, 0]
        for i in range(FREE_DIMENSION) for j in range(FREE_DIMENSION)
        for u in range(UNKNOWN_DIMENSION)
    )

    jets = []
    for column in range(85):
        gradient = [weights[column] * weight_direction[column][i] for i in range(FREE_DIMENSION)]
        hessian = [[weights[column] * (weight_direction[column][i] * weight_direction[column][j]
                    + sum((weight_u[column][u] * second[i][j][u, 0]
                           for u in range(UNKNOWN_DIMENSION)), arb(0)))
                    for j in range(FREE_DIMENSION)] for i in range(FREE_DIMENSION)]
        jets.append(Jet(weights[column], arb, gradient, hessian))
    parameters = json.loads(PARAMETERS.read_text(encoding="utf-8"))
    coefficients = [exact(fraction(x), arb, fmpq) for x in parameters["coefficient_rational_pairs"]]
    scales = [exact(fraction(x), arb, fmpq) for x in parameters["feature_scale_rational_pairs"]]
    direction = [exact(fraction(x), arb, fmpq) for x in parameters["shape_direction_rational_pairs"]]
    base_weights = [exact(fraction(x), arb, fmpq) for x in parameters["base_weight_rational_pairs"]]
    base_total = exact(fraction(parameters["base_total_weight_rational_pair"]), arb, fmpq)
    total = sum(jets, Jet(0, arb))
    quotient = (total - base_total) / (total + base_total)
    shape = sum((direction[c] * (jets[c] / total - base_weights[c] / base_total) ** 2
                 for c in range(85)), Jet(0, arb))
    correction = coefficients[0] * (scales[0] * quotient ** 2) ** 2 \
        + coefficients[1] * (scales[1] * shape) ** 2
    correction_metric = acb_mat(3, 3)
    for i in range(3):
        for j in range(3):
            correction_metric[i, j] = acb(
                (correction.hessian[i][j] + correction.hessian[i + 3][j + 3]) / 4,
                (correction.hessian[i][j + 3] - correction.hessian[i + 3][j]) / 4)
    correction_metric = (correction_metric + correction_metric.conjugate().transpose()) / 2
    correction_metric_max = maximum_finite_or_infinity(
        correction_metric[i, j]
        for i in range(3) for j in range(3)
    )

    tangent = [[acb(0) for _ in range(3)] for _ in range(85)]
    for local, index in enumerate(raw["free_indices"]):
        tangent[index][local] = acb(1)
        tangent[raw["eliminated"]][local] = acb(first_unknown[0, local], first_unknown[1, local])
    moment = arb_mat([[sum((raw["charges"][i][c] * raw["charges"][j][c] * weights[c]
                            for c in range(85)), arb(0)) for j in range(81)] for i in range(81)])
    vertical = acb_mat([[sum((raw["charges"][g][c] * weights[c] * tangent[c][j]
                             for c in range(85)), acb(0)) for j in range(3)] for g in range(81)])
    vertical_real_rhs = arb_mat([[vertical[i, j].real for j in range(3)] for i in range(81)])
    vertical_imag_rhs = arb_mat([[vertical[i, j].imag for j in range(3)] for i in range(81)])
    solved_real, moment_rho, moment_initial_norm, moment_tail = verified_solve(moment, vertical_real_rhs, [arb(1)] * 81, arb, arb_mat)
    solved_imag, _, _, _ = verified_solve(moment, vertical_imag_rhs, [arb(1)] * 81, arb, arb_mat)
    first_term = acb_mat([[sum((weights[c] * tangent[c][i].conjugate() * tangent[c][j]
                               for c in range(85)), acb(0)) for j in range(3)] for i in range(3)])
    moment_vertical = acb_mat([[acb(solved_real[i, j], solved_imag[i, j])
                                 for j in range(3)] for i in range(81)])
    baseline = first_term - vertical.conjugate().transpose() * moment_vertical
    baseline_max = maximum_finite_or_infinity(
        baseline[i, j]
        for i in range(3) for j in range(3)
    )
    corrected = (baseline + correction_metric + (baseline + correction_metric).conjugate().transpose()) / 2
    minors = []
    for dimension in (1, 2, 3):
        determinant = acb_mat([[corrected[i, j] for j in range(dimension)]
                               for i in range(dimension)]).det().real
        if not determinant.is_finite():
            minors.append({
                "dimension": dimension,
                "arb_enclosure": str(determinant),
                "finite_enclosure": False,
                "strictly_positive": False,
            })
            continue
        bounds = rational_bounds(determinant, arb, fmpq)
        minors.append({
            "dimension": dimension,
            "arb_enclosure": str(determinant),
            "finite_enclosure": True,
            "directed_rational_bounds": bounds,
            "lower_float": float(fraction(bounds["lower"])),
            "strictly_positive": fraction(bounds["lower"]) > 0,
        })
    primary = list(map(fraction, record["primary_leading_principal_minor_rational_lower_bounds"]))
    for independent, source in zip(minors, primary):
        if not independent["finite_enclosure"]:
            independent["primary_lower_inside_independent_enclosure"] = False
            continue
        bounds = independent["directed_rational_bounds"]
        independent["primary_lower_inside_independent_enclosure"] = (
            fraction(bounds["lower"]) <= source <= fraction(bounds["upper"])
        )
    return {
        "tube_record_index": record["tube_record_index"],
        "joint_neumann_rho_upper": float(joint_rho.upper()),
        "moment_neumann_rho_upper": float(moment_rho.upper()),
        "diagnostics": {
            "total_weight_lower_float": float(total.value.lower()),
            "total_weight_upper_float": float(total.value.upper()),
            "joint_solve_initial_scaled_rhs_norm": str(joint_initial_norm),
            "joint_solve_tail": str(joint_tail),
            "moment_solve_initial_rhs_norm": str(moment_initial_norm),
            "moment_solve_tail": str(moment_tail),
            "maximum_first_unknown_derivative_abs_upper": first_unknown_max,
            "full_verified_inverse_constructed": False,
            "nonfinite_first_unknown_derivative_count": len(nonfinite_first),
            "first_nonfinite_first_unknown_derivatives": nonfinite_first[:3],
            "maximum_second_unknown_derivative_abs_upper": second_unknown_max,
            "maximum_baseline_metric_entry_abs_upper": baseline_max,
            "maximum_correction_metric_entry_abs_upper": correction_metric_max,
        },
        "independent_leading_principal_minors": minors,
        "all_three_independent_minors_positive": all(x["strictly_positive"] for x in minors),
    }


def compute():
    from flint import acb, acb_mat, arb, arb_mat, ctx, fmpq
    ctx.dps = PRECISION_DIGITS
    inputs = json.loads(INPUTS.read_text(encoding="utf-8"))
    weakest_index = inputs["protocol"]["global_weakest_record_index"]
    selected = inputs["records"]
    if len(selected) != 8 or sum(item["tube_record_index"] == weakest_index for item in selected) != 1:
        raise ValueError("frozen stratified input selection changed")
    records = [verify_box(item, arb, acb, acb_mat, arb_mat, fmpq) for item in selected]
    passed = all(item["all_three_independent_minors_positive"] for item in records)
    return {
        "schema_version": 1,
        "inputs": {
            "generator_sha256": sha256(GENERATOR),
            "frozen_tube_inputs_sha256": sha256(INPUTS),
            "adapter_sha256": sha256(ADAPTER),
            "atlas_sha256": sha256(ATLAS),
            "exact_correction_parameters_sha256": sha256(PARAMETERS),
        },
        "protocol": {
            "selection_frozen_before_formula_reconstruction": True,
            "current_gate": "eight frozen stratified tube boxes including the global weakest",
            "imports_primary_metric_or_tube_generator": False,
            "common_primary_root_boxes_reused": True,
            "metric_rebuilt_on_serialized_krawczyk_image_not_wider_uniqueness_box": True,
            "six_real_free_coordinates_are_intervals": True,
            "quotient_metric_and_correction_hessian_reimplemented": True,
            "arb_precision_decimal_digits": PRECISION_DIGITS,
        },
        "records": records,
        "status": {
            "global_weakest_tube_box_independent_interval_metric_positive": passed,
            "all_eight_stratified_tube_boxes_independently_checked": len(records) == 8,
            "all_eight_stratified_tube_boxes_independent_interval_metrics_positive": len(records) == 8 and passed,
            "all_1084_tube_boxes_formula_independently_reimplemented": False,
            "root_box_existence_independently_reproved": False,
            "global_corrected_metric_positivity_certified": False,
        },
    }


def verify_payload(payload):
    if payload.get("schema_version") != 1 or len(payload.get("records", [])) != 8:
        return False
    expected_inputs = {
        "generator_sha256": sha256(GENERATOR),
        "frozen_tube_inputs_sha256": sha256(INPUTS),
        "adapter_sha256": sha256(ADAPTER),
        "atlas_sha256": sha256(ATLAS),
        "exact_correction_parameters_sha256": sha256(PARAMETERS),
    }
    if payload.get("inputs") != expected_inputs:
        return False
    frozen = json.loads(INPUTS.read_text(encoding="utf-8"))
    expected_indices = frozen["protocol"]["preselected_indices"]
    if [item.get("tube_record_index") for item in payload["records"]] != expected_indices:
        return False
    for record in payload["records"]:
        minors = record.get("independent_leading_principal_minors", [])
        if len(minors) != 3 or not record.get("all_three_independent_minors_positive"):
            return False
        for dimension, minor in enumerate(minors, 1):
            bounds = minor.get("directed_rational_bounds", {})
            if minor.get("dimension") != dimension or not minor.get("finite_enclosure"):
                return False
            if fraction(bounds["lower"]) <= 0 or fraction(bounds["lower"]) > fraction(bounds["upper"]):
                return False
            if not minor.get("strictly_positive"):
                return False
    protocol = payload.get("protocol", {})
    if protocol.get("imports_primary_metric_or_tube_generator"):
        return False
    if not protocol.get("six_real_free_coordinates_are_intervals"):
        return False
    status = payload.get("status", {})
    return bool(
        status.get("global_weakest_tube_box_independent_interval_metric_positive")
        and status.get("all_eight_stratified_tube_boxes_independently_checked")
        and status.get("all_eight_stratified_tube_boxes_independent_interval_metrics_positive")
        and not status.get("all_1084_tube_boxes_formula_independently_reimplemented")
        and not status.get("root_box_existence_independently_reproved")
        and not status.get("global_corrected_metric_positivity_certified")
    )


def main():
    payload = compute()
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
