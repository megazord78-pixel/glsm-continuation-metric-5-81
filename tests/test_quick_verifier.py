from __future__ import annotations

import copy
from fractions import Fraction
import json
import unittest

from glsmtube.endpoint_verify import ARTIFACT as ENDPOINT_ARTIFACT, verify_payload as verify_endpoint

from glsmtube.verify import (
    ATLAS, CertificateError, CORNER_TILING, EXACT_GEOMETRY, RAW_ADAPTER,
    ROBUST_CHAIN, ROOT_CERTIFICATE, TUBE_GRID, exact_integer, load_json,
    verify_exact_geometry, verify_package, verify_parent,
    verify_positive_minor_intervals, verify_robust_chain, verify_tiling,
)


class QuickVerifierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.roots = load_json(ROOT_CERTIFICATE)
        cls.tube = load_json(TUBE_GRID)
        cls.chain = load_json(ROBUST_CHAIN)
        cls.tiling = load_json(CORNER_TILING)

    def test_complete_release_certificate(self):
        result = verify_package()
        self.assertTrue(result["verified"])
        self.assertEqual(result["parent_root_boxes"], 1084)
        self.assertEqual(result["metric_child_boxes"], 34688)

    def test_component_margin_mutation_fails(self):
        record = copy.deepcopy(self.roots["records"][0])
        centre = record["unknown_center_rational_pairs"][2]
        radius = record["final_unknown_radius_rational_pairs"][2]
        c = Fraction(centre["numerator"], centre["denominator"])
        r = Fraction(radius["numerator"], radius["denominator"])
        bad = c - r
        record["dterm_krawczyk_image_directed_rational_bounds"][0]["lower"] = {
            "numerator": bad.numerator, "denominator": bad.denominator
        }
        with self.assertRaises(CertificateError):
            verify_parent(record, 0, self.tube["records"][0], self.chain)

    def test_derivative_origin_mutation_fails(self):
        record = copy.deepcopy(self.roots["records"][0])
        zero = {"numerator": 0, "denominator": 1}
        record["hypersurface_derivative_directed_rational_bounds"] = {
            "real": {"lower": zero, "upper": zero},
            "imag": {"lower": zero, "upper": zero},
        }
        with self.assertRaises(CertificateError):
            verify_parent(record, 0, self.tube["records"][0], self.chain)

    def test_corner_deletion_fails(self):
        changed = copy.deepcopy(self.tiling)
        changed["exact_tiling"]["corner_boxes"].pop()
        with self.assertRaises(CertificateError):
            verify_tiling(changed)

    def test_anchor_nesting_mutation_fails(self):
        changed_chain = copy.deepcopy(self.chain)
        record = self.roots["records"][0]
        patch = changed_chain["patches"][record["robust_patch_index"]]
        patch["krawczyk"]["residual_component_radius_rational_pairs"] = copy.deepcopy(
            record["final_unknown_radius_rational_pairs"]
        )
        with self.assertRaises(CertificateError):
            verify_parent(record, 0, self.tube["records"][0], changed_chain)

    def test_metric_lower_bound_mutation_fails(self):
        positive = {
            "lower": {"numerator": 1, "denominator": 10},
            "upper": {"numerator": 1, "denominator": 5},
        }
        changed = [copy.deepcopy(positive) for _ in range(3)]
        changed[2]["lower"] = {"numerator": 0, "denominator": 1}
        with self.assertRaises(CertificateError):
            verify_positive_minor_intervals(changed, 17)

    def test_nonintegral_geometry_is_never_rounded(self):
        with self.assertRaises(CertificateError):
            exact_integer(1.0000000000000002)

    def test_exact_geometry_registry_mutation_fails(self):
        registry = load_json(EXACT_GEOMETRY)
        adapter = load_json(RAW_ADAPTER)
        atlas = load_json(ATLAS)
        changed = copy.deepcopy(registry)
        changed["geometry"]["glsm_charges"][0][0] += 1
        with self.assertRaises(CertificateError):
            verify_exact_geometry(changed, adapter, atlas)

    def test_robust_patch_equality_mutation_fails(self):
        changed = copy.deepcopy(self.chain)
        patch = changed["patches"][0]["krawczyk"]
        patch["component_image_abs_upper_rational_pairs"][0] = copy.deepcopy(
            patch["residual_component_radius_rational_pairs"][0]
        )
        with self.assertRaises(CertificateError):
            verify_robust_chain(changed)

    def test_independent_endpoint_materiality_replay(self):
        result = verify_endpoint(load_json(ENDPOINT_ARTIFACT))
        self.assertTrue(result["nonempty_open_neighborhood_verified"])
        self.assertGreater(Fraction(result["largest_relative_correction_lower_exact"]), Fraction(1, 100))

    def test_endpoint_materiality_summary_mutation_fails(self):
        payload = load_json(ENDPOINT_ARTIFACT)
        changed = copy.deepcopy(payload)
        changed["replayed_evidence"]["materiality_witness_count"] += 1
        with self.assertRaises(ValueError):
            verify_endpoint(changed)


if __name__ == "__main__":
    unittest.main()
