from copy import deepcopy
import json
from pathlib import Path
import re
import unittest

from models.string_5_81_glsm_continuation_chain_independent_verifier import (
    ROOT, COVERS, BRIDGES, BRIDGE_3_COMPONENTS, ROBUST_CHAIN,
    BASE_ROOT, MATERIAL_BOX, sha256, verify_documents, verify_live,
)


def load_chain():
    covers = [json.loads((ROOT / "results" / name).read_text()) for name in COVERS]
    bridges = [json.loads((ROOT / "results" / name).read_text()) for name in BRIDGES]
    bridge3 = json.loads((ROOT / "results" / BRIDGE_3_COMPONENTS).read_text())
    robust = json.loads((ROOT / "results" / ROBUST_CHAIN).read_text())
    return covers, bridges, bridge3, robust


class ReleaseTests(unittest.TestCase):
    def test_live_publication_theorem(self):
        payload = verify_live()
        self.assertTrue(payload["status"]["publication_continuation_theorem_admitted"])
        self.assertTrue(payload["status"]["publication_transitive_provenance_release_admitted"])
        self.assertEqual(payload["proof_data"]["metric_segment_count"], 1003)
        self.assertEqual(payload["proof_data"]["patch_count"], 6)
        self.assertEqual(payload["proof_data"]["bridge_count"], 5)
        self.assertEqual(payload["proof_data"]["robust_root_chain"]["patch_count"], 7)
        self.assertEqual(payload["proof_data"]["robust_root_chain"]["bridge_count"], 7)
        self.assertFalse(payload["status"]["continuation_metric_formula_independently_reimplemented_on_all_1003_segments"])
        self.assertTrue(payload["status"]["continuation_metric_formula_independently_reimplemented_on_six_weakest_stratified_segments"])

    def test_deleted_segment_fails(self):
        covers, bridges, bridge3, robust = load_chain()
        covers[0]["segments"].pop()
        with self.assertRaises(ValueError):
            verify_documents(covers, bridges, bridge3, robust,
                             sha256(ROOT / "results" / BASE_ROOT),
                             sha256(ROOT / "results" / MATERIAL_BOX))

    def test_nonpositive_minor_fails(self):
        covers, bridges, bridge3, robust = load_chain()
        covers[0]["segments"][0]["leading_principal_minor_lower_bounds"][2]["real_lower"] = {
            "numerator": 0, "denominator": 1,
        }
        with self.assertRaises(ValueError):
            verify_documents(covers, bridges, bridge3, robust,
                             sha256(ROOT / "results" / BASE_ROOT),
                             sha256(ROOT / "results" / MATERIAL_BOX))

    def test_krawczyk_boundary_fails(self):
        covers, bridges, bridge3, robust = load_chain()
        record = bridges[0]["krawczyk"]["component_records"][0]
        record["image_abs_upper_rational"] = deepcopy(record["target_radius_rational"])
        with self.assertRaises(ValueError):
            verify_documents(covers, bridges, bridge3, robust,
                             sha256(ROOT / "results" / BASE_ROOT),
                             sha256(ROOT / "results" / MATERIAL_BOX))

    def test_manuscript_and_reviewer_docs_present(self):
        required = (
            "manuscript/main.tex",
            "manuscript/references.bib",
            "manuscript/main.bbl",
            "manuscript/main.pdf",
            "docs/TRUST_MODEL.md",
            "docs/CERTIFICATE_SCHEMA.md",
            "docs/PUBLICATION_CHECKLIST.md",
            "pyproject.toml",
            ".gitignore",
            ".gitattributes",
        )
        for relative in required:
            self.assertTrue((ROOT / relative).is_file(), relative)

    def test_manuscript_scope_is_fail_closed(self):
        manuscript = (ROOT / "manuscript" / "main.tex").read_text(encoding="utf-8")
        required_phrases = (
            "To the best of our knowledge",
            "not a global Ricci-flat metric",
            "or a complete string vacuum",
            "d\\in[1.49985,1.5]",
            "1003 exactly adjacent interval segments",
        )
        for phrase in required_phrases:
            self.assertIn(phrase, manuscript)
        self.assertNotIn("first corrected Calabi--Yau metric theorem", manuscript)

    def test_all_citations_have_bibliography_entries(self):
        manuscript = (ROOT / "manuscript" / "main.tex").read_text(encoding="utf-8")
        bibliography = (ROOT / "manuscript" / "references.bib").read_text(encoding="utf-8")
        cited = set()
        for group in re.findall(r"\\cite\{([^}]+)\}", manuscript):
            cited.update(item.strip() for item in group.split(","))
        entries = set(re.findall(r"@[A-Za-z]+\{([^,]+),", bibliography))
        self.assertEqual(set(), cited - entries)

    def test_complete_frozen_correction_is_packaged(self):
        path = ROOT / "data" / "frozen_correction_parameters.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        expected = {
            "coefficient": 2,
            "feature_scale": 2,
            "shape_direction": 85,
            "base_weight": 85,
        }
        for stem, count in expected.items():
            self.assertEqual(payload["dimensions"][f"{stem}_count"], count)
            self.assertEqual(len(payload[f"{stem}_rational_pairs"]), count)


if __name__ == "__main__":
    unittest.main()
