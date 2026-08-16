"""Independent interval reconstruction on the weakest segment of each patch.

The six segments are selected deterministically from the frozen covers by the
smallest serialized lower bound for the third leading principal minor.  This
module deliberately does not import any continuation metric generator or
adapter.  Certified root boxes are common inputs; the quotient metric,
implicit first/second derivatives, and correction Hessian are rebuilt here.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
ADAPTER = RESULTS / "2026-08-12_string_5_81_cymetric_toric_adapter.json"
ATLAS = RESULTS / "2026-08-12_string_5_81_d3_local_atlas_witness.json"
MATERIAL_BOX = RESULTS / "2026-08-13_string_5_81_glsm_depth15_material_interval_box.json"
SENSITIVITY = RESULTS / "2026-08-13_string_5_81_glsm_independent_joint_interval_derivative_verifier.json"
AFFINE_ROOT = RESULTS / "2026-08-13_string_5_81_glsm_affine_path_krawczyk_ladder.json"
ANISOTROPIC_ROOT = RESULTS / "2026-08-13_string_5_81_glsm_affine_path_root_anisotropic_refinement.json"
PARAMETERS = RESULTS / "2026-08-13_string_5_81_glsm_boundary_safe_exact_parameters.json"
RESULT = RESULTS / "2026-08-13_string_5_81_glsm_stratified_metric_formula_verifier.json"
PRECISION_DIGITS = 100
UNKNOWN_DIMENSION = 83
FREE_DIMENSION = 6
NEUMANN_ORDER = 50

PATCHES = [
    ("base_affine", "2026-08-13_string_5_81_glsm_affine_path_metric_connected_cover.json", None),
    ("recentered_1", "2026-08-13_string_5_81_glsm_recentered_patch_metric_cover.json",
     "2026-08-13_string_5_81_glsm_affine_path_recentered_root_patch.json"),
    ("recentered_2", "2026-08-13_string_5_81_glsm_recentered_patch_2_metric_cover.json",
     "2026-08-13_string_5_81_glsm_affine_path_recentered_root_patch_2.json"),
    ("local_tangent_1", "2026-08-13_string_5_81_glsm_local_tangent_metric_cover.json",
     "2026-08-13_string_5_81_glsm_local_tangent_root_patch.json"),
    ("local_tangent_2", "2026-08-13_string_5_81_glsm_local_tangent_2_metric_cover.json",
     "2026-08-13_string_5_81_glsm_local_tangent_root_patch_2.json"),
    ("local_tangent_3", "2026-08-13_string_5_81_glsm_local_tangent_3_metric_cover.json",
     "2026-08-13_string_5_81_glsm_local_tangent_root_patch_3.json"),
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fraction(value) -> Fraction:
    if isinstance(value, str):
        return Fraction(value)
    return Fraction(int(value["numerator"]), int(value["denominator"]))


def pair(value: Fraction) -> dict:
    return {"numerator": value.numerator, "denominator": value.denominator}


def exact(value: Fraction, arb, fmpq):
    return arb(fmpq(value.numerator, value.denominator))


def exact_integer(value) -> int:
    """Decode an integer-valued JSON scalar without rounding it.

    The upstream scientific input used JSON binary64 spellings such as
    ``1.0``.  Acceptance is therefore restricted to exactly integral values
    in the binary64 exact-integer range.  Nearby nonintegers fail closed.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("geometry entry is not a JSON integer/number")
    if isinstance(value, int):
        return value
    if not value.is_integer() or abs(value) > 2**53:
        raise ValueError("geometry entry is not an exactly represented integer")
    return int(value)


def rational_bounds(value, arb, fmpq) -> dict:
    import math
    lo = math.nextafter(float(value.lower()), -math.inf)
    while not arb(fmpq(*lo.as_integer_ratio())) < value:
        lo = math.nextafter(lo, -math.inf)
    hi = math.nextafter(float(value.upper()), math.inf)
    while not arb(fmpq(*hi.as_integer_ratio())) > value:
        hi = math.nextafter(hi, math.inf)
    lo_q = fmpq(*lo.as_integer_ratio())
    hi_q = fmpq(*hi.as_integer_ratio())
    return {
        "lower": {"numerator": int(lo_q.p), "denominator": int(lo_q.q)},
        "upper": {"numerator": int(hi_q.p), "denominator": int(hi_q.q)},
    }


def parse_acb(text, arb, acb):
    match = re.fullmatch(r"\(([-+0-9.eE]+)\s*([+-])\s*([0-9.eE+-]+)j\)", text.strip())
    if match is None:
        raise ValueError(f"invalid complex coefficient {text!r}")
    imaginary = arb(match.group(3))
    if match.group(2) == "-":
        imaginary = -imaginary
    return acb(arb(match.group(1)), imaginary)


class Jet:
    def __init__(self, value, arb_type, gradient=None, hessian=None):
        self.arb = arb_type
        self.value = value
        self.gradient = gradient or [arb_type(0) for _ in range(FREE_DIMENSION)]
        self.hessian = hessian or [[arb_type(0) for _ in range(FREE_DIMENSION)]
                                   for _ in range(FREE_DIMENSION)]

    def _jet(self, other):
        return other if isinstance(other, Jet) else Jet(other, self.arb)

    def __add__(self, other):
        other = self._jet(other)
        return Jet(self.value + other.value, self.arb,
                   [a + b for a, b in zip(self.gradient, other.gradient)],
                   [[self.hessian[i][j] + other.hessian[i][j]
                     for j in range(FREE_DIMENSION)] for i in range(FREE_DIMENSION)])

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.value, self.arb, [-x for x in self.gradient],
                   [[-x for x in row] for row in self.hessian])

    def __sub__(self, other):
        return self + (-self._jet(other))

    def __rsub__(self, other):
        return self._jet(other) - self

    def __mul__(self, other):
        other = self._jet(other)
        gradient = [self.gradient[i] * other.value + self.value * other.gradient[i]
                    for i in range(FREE_DIMENSION)]
        hessian = [[
            self.hessian[i][j] * other.value
            + self.gradient[i] * other.gradient[j]
            + self.gradient[j] * other.gradient[i]
            + self.value * other.hessian[i][j]
            for j in range(FREE_DIMENSION)] for i in range(FREE_DIMENSION)]
        return Jet(self.value * other.value, self.arb, gradient, hessian)

    __rmul__ = __mul__

    def reciprocal(self):
        inv = 1 / self.value
        gradient = [-self.gradient[i] * inv ** 2 for i in range(FREE_DIMENSION)]
        hessian = [[
            2 * self.gradient[i] * self.gradient[j] * inv ** 3
            - self.hessian[i][j] * inv ** 2
            for j in range(FREE_DIMENSION)] for i in range(FREE_DIMENSION)]
        return Jet(inv, self.arb, gradient, hessian)

    def __truediv__(self, other):
        return self * self._jet(other).reciprocal()

    def __rtruediv__(self, other):
        return self._jet(other) / self

    def __pow__(self, exponent):
        if exponent < 0:
            return (self.reciprocal()) ** (-exponent)
        result = Jet(1, self.arb)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent >>= 1
        return result


def verified_inverse(matrix, component_radii, arb, arb_mat):
    midpoint = arb_mat([[matrix[i, j].mid() for j in range(matrix.ncols())]
                        for i in range(matrix.nrows())])
    preconditioner = midpoint.inv()
    identity = arb_mat(matrix.nrows(), matrix.ncols())
    for i in range(matrix.nrows()):
        identity[i, i] = 1
    remainder = identity - preconditioner * matrix
    scaled = arb_mat([[
        remainder[i, j] * component_radii[j] / component_radii[i]
        for j in range(matrix.ncols())] for i in range(matrix.nrows())])
    rho = max(sum((abs(scaled[i, j]) for j in range(matrix.ncols())), arb(0))
              for i in range(matrix.nrows()))
    if not rho < 1:
        raise ValueError(f"independent scaled Neumann norm is not below one: {rho}")
    right = arb_mat([[preconditioner[i, j] / component_radii[i]
                      for j in range(matrix.ncols())] for i in range(matrix.nrows())])
    right_norm = max(sum((abs(right[i, j]) for j in range(matrix.ncols())), arb(0))
                     for i in range(matrix.nrows()))
    series = identity
    term = identity
    for _ in range(NEUMANN_ORDER):
        term = term * scaled
        series = series + term
    approximate = series * right
    tail = rho ** (NEUMANN_ORDER + 1) / (1 - rho) * right_norm
    inverse = arb_mat([[
        component_radii[i] * approximate[i, j] + arb(0, component_radii[i] * tail)
        for j in range(matrix.ncols())] for i in range(matrix.nrows())])
    return inverse, rho


def raw_inputs(arb, acb, fmpq):
    raw = json.loads(ADAPTER.read_text(encoding="utf-8"))["toric_sampler_core"]
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    center = json.loads(MATERIAL_BOX.read_text(encoding="utf-8"))["exact_center"]
    sensitivity = json.loads(SENSITIVITY.read_text(encoding="utf-8"))
    selected = [int(value) - 1 for value in atlas["fan"]["selected_cone_divisor_indices"]]
    y0_pairs = center["eliminated_center_rational_pairs"] + center["dterm_center_rational_pairs"]
    slope_pairs = []
    for bounds in sensitivity["path_root_sensitivity"]["depth_derivative_directed_rational_bounds"]:
        slope_pairs.append(pair((fraction(bounds["lower"]) + fraction(bounds["upper"])) / 2))
    return {
        "charges": [[exact_integer(x) for x in row] for row in raw["glsm_charges"]],
        "exponents": [[exact_integer(x) for x in row] for row in raw["exp_aK"]],
        "coefficients": [parse_acb(x, arb, acb) for x in raw["coeff_aK_text"]],
        "free_indices": [selected[0], selected[1], selected[3]],
        "free_index": selected[0],
        "eliminated": selected[2],
        "base_intercept_pairs": y0_pairs,
        "base_slope_pairs": slope_pairs,
    }


def selected_segment(cover):
    return min(cover["segments"],
               key=lambda item: fraction(item["leading_principal_minor_lower_bounds"][2]["real_lower"]))


def t_interval(segment):
    if "physical_t_lower" in segment:
        return fraction(segment["physical_t_lower"]), fraction(segment["physical_t_upper"])
    return (fraction(segment["depth_lower"]) - Fraction(3, 2),
            fraction(segment["depth_upper"]) - Fraction(3, 2))


def root_configuration(name, segment, root_name):
    lower, upper = t_interval(segment)
    if name == "base_affine":
        ladder = json.loads(AFFINE_ROOT.read_text(encoding="utf-8"))
        anisotropic = json.loads(ANISOTROPIC_ROOT.read_text(encoding="utf-8"))
        required = max(abs(lower), abs(upper))
        choices = [record for record in ladder["records"]
                   if record["certified"] and Fraction(record["path_radius"]) >= required]
        target = min(choices, key=lambda record: Fraction(record["path_radius"]))
        selected = anisotropic["records"].get(target["path_radius"], target["selected"])
        radii = selected.get("residual_component_radius_rational_pairs")
        if radii is None:
            radii = [selected["residual_component_radius_rational_pair"]] * UNKNOWN_DIMENSION
        return None, None, radii, [AFFINE_ROOT, ANISOTROPIC_ROOT]
    path = RESULTS / root_name
    root = json.loads(path.read_text(encoding="utf-8"))
    if name.startswith("recentered"):
        center = root["center"]["residual_center_rational_pairs"]
        radii = root["krawczyk"]["residual_component_radius_rational_pairs"]
        return "base", center, radii, [path]
    record = next(item for item in root["records"]
                  if fraction(item["u_radius"]) == Fraction("3e-5"))
    if not record["certified"]:
        raise ValueError("selected local root record is not certified")
    return root["local_affine_coordinates"]["intercept_rational_pairs"], None, \
        record["krawczyk"]["residual_component_radius_rational_pairs"], [path]


def verify_segment(name, cover_path, root_name, arb, acb, acb_mat, arb_mat, fmpq):
    cover = json.loads(cover_path.read_text(encoding="utf-8"))
    segment = selected_segment(cover)
    lower, upper = t_interval(segment)
    raw = raw_inputs(arb, acb, fmpq)
    intercept_spec, offset_pairs, radius_pairs, root_paths = root_configuration(name, segment, root_name)
    if intercept_spec is None or intercept_spec == "base":
        intercept_pairs = raw["base_intercept_pairs"]
        slope_pairs = raw["base_slope_pairs"]
    else:
        intercept_pairs = intercept_spec
        root = json.loads(root_paths[0].read_text(encoding="utf-8"))
        slope_pairs = root["local_affine_coordinates"]["slope_rational_pairs"]
    intercept = [exact(fraction(x), arb, fmpq) for x in intercept_pairs]
    if offset_pairs is not None:
        intercept = [x + exact(fraction(y), arb, fmpq) for x, y in zip(intercept, offset_pairs)]
    slopes = [exact(fraction(x), arb, fmpq) for x in slope_pairs]
    radii = [exact(fraction(x), arb, fmpq) for x in radius_pairs]
    t = arb(exact((lower + upper) / 2, arb, fmpq), exact((upper - lower) / 2, arb, fmpq))
    unknown = [intercept[i] + slopes[i] * t + arb(0, radii[i]) for i in range(UNKNOWN_DIMENSION)]

    logs = [acb(0) for _ in range(85)]
    logs[raw["free_index"]] = acb(arb("-1.5") - t)
    logs[raw["eliminated"]] = acb(unknown[0], unknown[1])
    terms = []
    for coefficient, exponent in zip(raw["coefficients"], raw["exponents"]):
        argument = sum((exponent[i] * logs[i] for i in range(85)), acb(0))
        terms.append(coefficient * argument.exp())
    weights = []
    for column in range(85):
        log_abs = (2 * unknown[0] if column == raw["eliminated"] else arb(0))
        if column == raw["free_index"]:
            log_abs = 2 * (arb("-1.5") - t)
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
    inverse, joint_rho = verified_inverse(jacobian, radii, arb, arb_mat)
    first_unknown = -(inverse * free_derivative)
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
            second[i][j] = second[j][i] = -(inverse * rhs)

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
    q = (total - base_total) / (total + base_total)
    shape = sum((direction[c] * (jets[c] / total - base_weights[c] / base_total) ** 2
                 for c in range(85)), Jet(0, arb))
    correction = coefficients[0] * (scales[0] * q ** 2) ** 2 \
        + coefficients[1] * (scales[1] * shape) ** 2
    correction_metric = acb_mat(3, 3)
    for i in range(3):
        for j in range(3):
            correction_metric[i, j] = acb(
                (correction.hessian[i][j] + correction.hessian[i + 3][j + 3]) / 4,
                (correction.hessian[i][j + 3] - correction.hessian[i + 3][j]) / 4)
    correction_metric = (correction_metric + correction_metric.conjugate().transpose()) / 2

    tangent = [[acb(0) for _ in range(3)] for _ in range(85)]
    for local, index in enumerate(raw["free_indices"]):
        tangent[index][local] = acb(1)
        tangent[raw["eliminated"]][local] = acb(first_unknown[0, local], first_unknown[1, local])
    moment = arb_mat([[sum((raw["charges"][i][c] * raw["charges"][j][c] * weights[c]
                            for c in range(85)), arb(0)) for j in range(81)] for i in range(81)])
    moment_inverse, moment_rho = verified_inverse(moment, [arb(1)] * 81, arb, arb_mat)
    first_term = acb_mat([[sum((weights[c] * tangent[c][i].conjugate() * tangent[c][j]
                               for c in range(85)), acb(0)) for j in range(3)] for i in range(3)])
    vertical = acb_mat([[sum((raw["charges"][g][c] * weights[c] * tangent[c][j]
                             for c in range(85)), acb(0)) for j in range(3)] for g in range(81)])
    moment_complex = acb_mat([[acb(moment_inverse[i, j]) for j in range(81)] for i in range(81)])
    baseline = first_term - vertical.conjugate().transpose() * moment_complex * vertical
    corrected = (baseline + correction_metric + (baseline + correction_metric).conjugate().transpose()) / 2
    minors = []
    for dimension in (1, 2, 3):
        determinant = acb_mat([[corrected[i, j] for j in range(dimension)]
                               for i in range(dimension)]).det().real
        bounds = rational_bounds(determinant, arb, fmpq)
        minors.append({"dimension": dimension, "directed_rational_bounds": bounds,
                       "lower_float": float(fraction(bounds["lower"])),
                       "strictly_positive": fraction(bounds["lower"]) > 0})
    primary = segment["leading_principal_minor_lower_bounds"]
    for independent, source in zip(minors, primary):
        source_lower = fraction(source["real_lower"])
        enclosure = independent["directed_rational_bounds"]
        independent["primary_lower_inside_independent_enclosure"] = (
            fraction(enclosure["lower"]) <= source_lower <= fraction(enclosure["upper"])
        )
        independent["relative_lower_difference_from_primary"] = float(
            abs(fraction(enclosure["lower"]) - source_lower) / source_lower
        )
    return {
        "patch": name,
        "cover_sha256": sha256(cover_path),
        "root_input_sha256": {path.name: sha256(path) for path in root_paths},
        "selection": {
            "rule": "minimum exact serialized lower bound of third leading principal minor",
            "segment_index": segment["index"],
            "physical_t_interval": [pair(lower), pair(upper)],
            "primary_third_minor_lower": primary[2]["real_lower"],
            "primary_third_minor_lower_float": primary[2]["real_lower_float"],
        },
        "inverse_bounds": {"joint_neumann_rho_upper": float(joint_rho.upper()),
                           "moment_neumann_rho_upper": float(moment_rho.upper())},
        "independent_leading_principal_minors": minors,
        "all_three_independent_minors_positive": all(x["strictly_positive"] for x in minors),
    }


def compute():
    from flint import acb, acb_mat, arb, arb_mat, ctx, fmpq
    ctx.dps = PRECISION_DIGITS
    records = []
    for name, cover_name, root_name in PATCHES:
        records.append(verify_segment(name, RESULTS / cover_name, root_name,
                                      arb, acb, acb_mat, arb_mat, fmpq))
    passed = all(record["all_three_independent_minors_positive"] for record in records)
    all_primary_lowers_enclosed = all(
        minor["primary_lower_inside_independent_enclosure"]
        for record in records for minor in record["independent_leading_principal_minors"]
    )
    return {
        "schema_version": 1,
        "protocol": {
            "selection_frozen_before_independent_reconstruction": True,
            "selection_rule": "one exact weakest-third-minor segment from each of six continuation patches",
            "arb_precision_decimal_digits": PRECISION_DIGITS,
            "imports_primary_metric_generator_or_adapter": False,
            "common_certified_root_boxes_reused": True,
            "quotient_metric_and_correction_formula_reimplemented": True,
        },
        "inputs": {path.name: sha256(path) for path in
                   (ADAPTER, ATLAS, MATERIAL_BOX, SENSITIVITY, PARAMETERS)},
        "records": records,
        "status": {
            "six_deterministically_stratified_segments_checked": len(records) == 6,
            "all_six_independent_interval_metrics_positive": passed,
            "stratified_metric_formula_replication_passed": passed,
            "all_18_primary_minor_lowers_inside_independent_enclosures": all_primary_lowers_enclosed,
            "all_1003_metric_segments_formula_independently_reimplemented": False,
            "root_box_existence_independently_reproved": False,
            "global_corrected_metric_positivity_certified": False,
        },
        "scope": [
            "This independently reconstructs the formula on six preselected weak segments.",
            "It does not independently reprove the shared rational root boxes.",
            "It does not extend formula independence to the other 997 segments or to the full manifold.",
        ],
    }


def verify_payload(payload: dict) -> bool:
    """Fast final verifier for selection, provenance, rational signs, and scope."""
    if payload.get("schema_version") != 1 or len(payload.get("records", [])) != 6:
        return False
    expected_names = [item[0] for item in PATCHES]
    if [item.get("patch") for item in payload["records"]] != expected_names:
        return False
    for (name, cover_name, root_name), record in zip(PATCHES, payload["records"]):
        cover_path = RESULTS / cover_name
        if record.get("cover_sha256") != sha256(cover_path):
            return False
        chosen = selected_segment(json.loads(cover_path.read_text(encoding="utf-8")))
        if record.get("selection", {}).get("segment_index") != chosen["index"]:
            return False
        if len(record.get("independent_leading_principal_minors", [])) != 3:
            return False
        for minor in record["independent_leading_principal_minors"]:
            bounds = minor.get("directed_rational_bounds", {})
            if fraction(bounds["lower"]) <= 0 or fraction(bounds["lower"]) > fraction(bounds["upper"]):
                return False
            if not minor.get("primary_lower_inside_independent_enclosure"):
                return False
        if not record.get("all_three_independent_minors_positive"):
            return False
    inputs = payload.get("inputs", {})
    for path in (ADAPTER, ATLAS, MATERIAL_BOX, SENSITIVITY, PARAMETERS):
        if inputs.get(path.name) != sha256(path):
            return False
    status = payload.get("status", {})
    return bool(
        status.get("stratified_metric_formula_replication_passed")
        and status.get("all_18_primary_minor_lowers_inside_independent_enclosures")
        and not status.get("all_1003_metric_segments_formula_independently_reimplemented")
        and not status.get("root_box_existence_independently_reproved")
        and not status.get("global_corrected_metric_positivity_certified")
    )


def main():
    payload = compute()
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"records": [{"patch": x["patch"], "segment": x["selection"]["segment_index"],
                                    "minors": x["independent_leading_principal_minors"]}
                                   for x in payload["records"]], "status": payload["status"]}, indent=2))
    raise SystemExit(0 if verify_payload(payload) else 1)


if __name__ == "__main__":
    main()
