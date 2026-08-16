"""Cached, append-only, resumable replay of 32 corner boxes at all 1084 path centres."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
from itertools import product
import json
import os
from pathlib import Path

from regeneration.models.string_5_81_glsm_cached_interval_engine import install_worker_caches
from regeneration.models.string_5_81_glsm_transverse_tube_full_independent_replay import (
    ADAPTER, ATLAS, CHAIN, INITIAL_UNKNOWN_RADIUS, TUBE,
    initial_record, pair, refined_verify_record,
)


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = Path(__file__).resolve()
ENGINE = ROOT / "models/string_5_81_glsm_cached_interval_engine.py"
WORK = ROOT / "work"
TILING = ROOT / "results/2026-08-14_string_5_81_glsm_full_path_cube_corner_tiling_protocol.json"
EQUIVALENCE = ROOT / "results/2026-08-14_string_5_81_glsm_cached_engine_equivalence_benchmark.json"
FACTOR2 = ROOT / "results/2026-08-14_string_5_81_glsm_stratified_factor2_weak_vertex_screen.json"
RESULT = WORK / "metric_regenerated.json"
CHECKPOINT_META = WORK / "metric_checkpoint_meta.json"
CHECKPOINT_JSONL = WORK / "metric_checkpoint.jsonl"
LONGITUDINAL_RADIUS = Fraction(7, 100_000_000)
TRANSVERSE_RADIUS = Fraction(1, 200_000_000)
WORKERS = max(1, int(os.environ.get("GLSM_FULL_CUBE_WORKERS", "8")))

_TUBE_PAYLOAD = None
_CHAIN_PAYLOAD = None
_ADAPTER_PAYLOAD = None
_ATLAS_PAYLOAD = None
_VERIFY_BOX = None


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def q(value):
    return Fraction(value["numerator"], value["denominator"])


def input_hashes():
    return {"generator_sha256": sha256(GENERATOR), "cached_engine_sha256": sha256(ENGINE),
            "primary_tube_grid_sha256": sha256(TUBE), "robust_krawczyk_chain_sha256": sha256(CHAIN),
            "corner_tiling_protocol_sha256": sha256(TILING), "cached_equivalence_benchmark_sha256": sha256(EQUIVALENCE),
            "factor2_screen_sha256": sha256(FACTOR2)}


def initialize_worker():
    global _TUBE_PAYLOAD, _CHAIN_PAYLOAD, _ADAPTER_PAYLOAD, _ATLAS_PAYLOAD, _VERIFY_BOX
    install_worker_caches()
    _TUBE_PAYLOAD = json.loads(TUBE.read_text(encoding="utf-8"))
    _CHAIN_PAYLOAD = json.loads(CHAIN.read_text(encoding="utf-8"))
    _ADAPTER_PAYLOAD = json.loads(ADAPTER.read_text(encoding="utf-8"))["toric_sampler_core"]
    _ATLAS_PAYLOAD = json.loads(ATLAS.read_text(encoding="utf-8"))
    from regeneration.models.string_5_81_glsm_transverse_tube_independent_formula_verifier import verify_box
    _VERIFY_BOX = verify_box


def task(spec):
    from flint import acb, acb_mat, arb, arb_mat, ctx, fmpq
    ctx.dps = 100
    sequence, path_index, slots = spec
    record = initial_record(_TUBE_PAYLOAD["records"][path_index], _CHAIN_PAYLOAD["patches"])
    centers = [q(item) for item in record["free_center_rational_pairs"]]
    for coordinate, slot in enumerate(slots, start=1):
        centers[coordinate] += slot * TRANSVERSE_RADIUS
    record["free_center_rational_pairs"] = [pair(value) for value in centers]
    record["free_radius_rational_pairs"] = [pair(LONGITUDINAL_RADIUS)] + [pair(TRANSVERSE_RADIUS)] * 5
    record["unknown_radius_rational_pairs"] = [pair(INITIAL_UNKNOWN_RADIUS)] * 83
    root = refined_verify_record(record, _ADAPTER_PAYLOAD, _ATLAS_PAYLOAD, arb, acb, arb_mat, fmpq)
    grid = _TUBE_PAYLOAD["records"][path_index]
    base = {"sequence": sequence, "path_index": path_index, "offset_slots": list(slots),
            "t_center": grid["t_center"], "depth_center": grid["depth_center"],
            "root_history": root["independent_refinement_history"],
            "root_pass": root["root_existence_independently_reproved"]}
    if not base["root_pass"]:
        return sequence, {**base, "metric_pass": False, "full_pass": False}
    metric_record = dict(record)
    metric_record["unknown_radius_rational_pairs"] = root["final_unknown_radius_rational_pairs"]
    metric_record["primary_hypersurface_newton_image_directed_rational_bounds"] = root["hypersurface"]["newton_image_directed_rational_bounds"]
    metric_record["primary_dterm_krawczyk_image_directed_rational_bounds"] = root["dterm"]["image_directed_rational_bounds"]
    metric = _VERIFY_BOX(metric_record, arb, acb, acb_mat, arb_mat, fmpq)
    metric_pass = metric["all_three_independent_minors_positive"]
    return sequence, {**base, "root_q": root["dterm"]["maximum_image_radius_ratio"],
                      "minor_directed_rational_bounds": [item["directed_rational_bounds"] for item in metric["independent_leading_principal_minors"]],
                      "metric_pass": metric_pass, "full_pass": metric_pass}


def load_checkpoint(expected, total):
    if CHECKPOINT_META.exists():
        meta = json.loads(CHECKPOINT_META.read_text(encoding="utf-8"))
        if meta.get("inputs") != expected or meta.get("total_record_count") != total:
            raise ValueError("checkpoint provenance/protocol mismatch")
    elif CHECKPOINT_JSONL.exists():
        raise ValueError("checkpoint JSONL exists without provenance meta")
    else:
        CHECKPOINT_META.write_text(json.dumps({"schema_version": 1, "inputs": expected,
                                              "total_record_count": total, "format": "one canonical JSON record per line"}, indent=2) + "\n",
                                   encoding="utf-8", newline="\n")
    records = {}
    truncated_final_line_ignored = False
    if CHECKPOINT_JSONL.exists():
        lines = CHECKPOINT_JSONL.read_text(encoding="utf-8").splitlines()
        for line_index, line in enumerate(lines):
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                if line_index != len(lines) - 1:
                    raise ValueError(f"malformed non-final checkpoint line {line_index}")
                truncated_final_line_ignored = True
                continue
            sequence = record.get("sequence")
            if not isinstance(sequence, int) or not 0 <= sequence < total or sequence in records:
                raise ValueError(f"invalid/duplicate checkpoint sequence {sequence}")
            records[sequence] = record
    return records, truncated_final_line_ignored


def compute():
    tiling = json.loads(TILING.read_text(encoding="utf-8"))
    equivalence = json.loads(EQUIVALENCE.read_text(encoding="utf-8"))
    factor2 = json.loads(FACTOR2.read_text(encoding="utf-8"))
    if not tiling["status"]["exact_2pow5_corner_tiling_protocol_frozen"]:
        raise ValueError("corner tiling not frozen")
    if not equivalence["status"]["cached_engine_admitted_for_full_path_replay"]:
        raise ValueError("cached engine not admitted")
    if factor2["status"]["factor2_all_fifteen_strata_weak_vertices_certified"]:
        raise ValueError("factor-two screen no longer justifies base radius")
    tube = json.loads(TUBE.read_text(encoding="utf-8"))
    vertices = list(product((-1, 1), repeat=5))
    specs = [(path_index * 32 + vertex_index, path_index, slots)
             for path_index in range(len(tube["records"])) for vertex_index, slots in enumerate(vertices)]
    expected = input_hashes()
    records, truncated = load_checkpoint(expected, len(specs))
    missing = [spec for spec in specs if spec[0] not in records]
    if missing:
        with CHECKPOINT_JSONL.open("a", encoding="utf-8", newline="\n", buffering=1) as checkpoint:
            with ProcessPoolExecutor(max_workers=min(WORKERS, len(missing)), initializer=initialize_worker) as executor:
                for sequence, record in executor.map(task, missing, chunksize=1):
                    checkpoint.write(json.dumps(record, separators=(",", ":")) + "\n")
                    records[sequence] = record
                    if not record["full_pass"]:
                        raise RuntimeError(f"full-path corner box failed at sequence {sequence}")
                    if len(records) % 100 == 0:
                        print(json.dumps({"completed": len(records), "total": len(specs)}), flush=True)
    ordered = [records[index] for index in range(len(specs))]
    minimum = min(q(item["minor_directed_rational_bounds"][2]["lower"]) for item in ordered)
    return {
        "schema_version": 1, "inputs": expected,
        "protocol": {"path_center_count": len(tube["records"]), "corner_count_per_path_center": 32,
                     "record_count": len(ordered), "corner_alphabet": [-1, 1],
                     "longitudinal_radius": pair(LONGITUDINAL_RADIUS), "transverse_radius": pair(TRANSVERSE_RADIUS),
                     "worker_count": WORKERS, "cached_engine_enabled": True,
                     "checkpoint_meta": str(CHECKPOINT_META.relative_to(ROOT)).replace("\\", "/"),
                     "checkpoint_jsonl": str(CHECKPOINT_JSONL.relative_to(ROOT)).replace("\\", "/"),
                     "truncated_final_checkpoint_line_ignored": truncated},
        "records": ordered,
        "summary": {"root_pass_count": sum(item["root_pass"] for item in ordered),
                    "metric_pass_count": sum(item["metric_pass"] for item in ordered),
                    "minimum_third_minor_lower_exact": str(minimum),
                    "maximum_root_q": max(item["root_q"] for item in ordered)},
        "status": {"all_34688_corner_boxes_independently_certified": all(item["full_pass"] for item in ordered),
                   "full_path_transverse_cube_radius_1e_8_certified": all(item["full_pass"] for item in ordered),
                   "global_corrected_metric_positivity_on_cy_certified": False,
                   "publication_release_admitted": False},
    }


def main():
    payload = compute()
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"protocol": payload["protocol"], "summary": payload["summary"], "status": payload["status"]}, indent=2))


if __name__ == "__main__":
    main()
