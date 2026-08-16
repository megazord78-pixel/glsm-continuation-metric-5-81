# Interval-certified corrected GLSM metric on CY `(5,81)`

This standalone `v2.0.0` research package consolidates the earlier central-path
calculation, its independently reconstructed endpoint, and the product-tube
extension into one computer-assisted theorem hierarchy:

> On the disclosed GLSM chart and for the exactly specified gauge-invariant
> Kähler ansatz, there is one locally unique parameterized implicit root branch,
> and the corrected Hermitian metric is positive on
>
> `t in [-3/20000, 0]` and `u in [-1/100000000, 1/100000000]^5`.

The package separates two verification levels:

- `python verify_release.py --quick` authenticates the immutable release,
  verifies all exact rational root margins, the complete product tiling,
  longitudinal branch gluing, and all stored directed metric-minor intervals;
- `python verify_release.py --full` additionally regenerates the Arb root and
  metric enclosures from the frozen formula inputs. This mode requires
  `python-flint==0.9.0` / FLINT `3.6.0` and is intentionally much slower.

The quick verifier is a last-mile exact checker. It does not claim to regenerate
the Arb formula enclosures.

The full Arb regeneration has completed all 1,084 root parents and all 34,688
metric children, with exact equality of every proof-relevant field to the
archived release.  The quick verifier additionally replays an independent
83-real endpoint root and metric reconstruction and a fixed 31-direction
materiality test.  Two directions have a certified relative correction above
one percent; the largest lower bound is `0.0183143620836...`.

## Certified domain

The longitudinal path is

`d in [1.49985, 1.5]`, equivalently `t=d-1.5 in [-3/20000,0]`.

Five real transverse coordinates satisfy

`u_i in [-1e-8,1e-8]`.

There are 1,084 parent root boxes and 34,688 metric child boxes. The 32
children at each path centre exactly tile the five-dimensional transverse cube.

## Scope boundary

This package does **not** certify global positivity on the compact
Calabi--Yau, Ricci flatness, approximation error to the canonical Ricci-flat
metric, a complete physical alpha-prime correction, or a physical string
vacuum.

## Layout

- `data/` — compressed proof records and exact frozen inputs;
- `data/endpoint_validation/` — independent endpoint and materiality evidence;
- `src/glsmtube/` — dependency-free quick verifier;
- `regeneration/` — optional Arb formula regeneration;
- `tests/` — positive and fail-closed mutation tests;
- `docs/` — certificate schema, gluing lemma, trust model and reproducibility;
- `manuscript/` — article source and compiled PDF;
- `MANIFEST.sha256.json` — immutable byte-level release binding.

## Citation

Citation metadata are provided in `CITATION.cff`. Code is BSD-3-Clause. The
archival DOI for version 2.0.0 will be inserted only after the immutable release
has been deposited.
