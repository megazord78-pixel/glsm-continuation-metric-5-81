# Interval-certified corrected GLSM metric on the `(5,81)` model

This is the standalone `v1.0.1` research release accompanying:

> Evgeniy Agafonov, *An Interval-Certified Positivity Theorem for a
> Corrected GLSM Metric along a Calabi--Yau Continuation Path* (2026).

The central computer-assisted theorem is

\[
d\in[1.49985,1.5],\qquad G_{\mathrm{corr}}(d)\succ0.
\]

Here `corrected` means the exact frozen gauge-invariant `ddbar` ansatz stated
in the manuscript.  It does not mean that a complete physical alpha-prime
correction has been derived.

## Evidence at a glance

- one 83-real implicit branch covered by 7 Krawczyk patches and 7 bridges;
- every robust componentwise Krawczyk ratio is strictly below `0.85`;
- 6 metric patches and 1003 exactly adjacent interval segments;
- three strictly positive leading-principal-minor lower bounds per segment;
- 18 independently reconstructed positive minors on the weakest segment of
  each metric patch;
- an independent endpoint root and exact-rational determinant replay;
- a fixed-direction relative correction lower bound above `1.83%`;
- zero unresolved or ambiguous SHA-256 edges in the proof provenance graph.

## Verify the release

Python 3.11 or later is the only requirement:

```powershell
python verify_release.py
```

The command authenticates every release file, runs the dependency-free
central checker, and executes mutation tests.  It does not regenerate the
heavy Arb enclosures.

An installed package also exposes:

```powershell
verify-glsm-metric-certificate
```

## Repository structure

- `manuscript/main.tex` — complete article;
- `manuscript/references.bib` — bibliography;
- `results/` — exact interval and provenance records;
- `data/frozen_correction_parameters.json` — the complete exact-rational
  frozen correction ansatz (2 coefficients, 2 feature scales, 85 direction
  entries, and 85 reference weights);
- `models/string_5_81_glsm_continuation_chain_independent_verifier.py` —
  stdlib-only proof replay;
- `tests/` — positive and negative mutation tests;
- `docs/TRUST_MODEL.md` — trusted base and evidence boundary;
- `docs/CERTIFICATE_SCHEMA.md` — proof-record format;
- `docs/PUBLICATION_CHECKLIST.md` — public-release and submission sequence;
- `MANIFEST.sha256.json` — complete byte-level release binding.

## Scope boundary

The theorem concerns a one-dimensional trajectory in a disclosed toric
hypersurface chart.  It does not prove positivity on an open tube, positivity
on the full compact hypersurface, global Ricci flatness, a global
Monge--Ampere error bound, or a complete physical string vacuum.

The 1003-segment layer replays serialized Arb bounds exactly.  A second
formula implementation covers one deterministically selected weak segment
from each of the six patches, not all 1003 segments, and reuses the certified
root boxes.  These independence boundaries are part of the scientific claim.

## Citation and license

Citation metadata are in `CITATION.cff`.  Code and verification material are
distributed under the BSD 3-Clause License.  Release `v1.0.1` is an
archive-integration update: it does not alter the theorem, proof records, or
numerical bounds.  The Zenodo DOI will be added after ingestion completes.
