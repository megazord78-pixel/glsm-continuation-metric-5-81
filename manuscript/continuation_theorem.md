# A computer-assisted positivity theorem along a GLSM continuation path

Evgeniy Agafonov — Independent Researcher, Sochi, Russia  
ORCID: 0009-0005-9059-5291

## Result

For the disclosed toric hypersurface chart and frozen correction ansatz, the
certified solution branch exists and is locally unique along
`d in [1.49985, 1.5]`. The corrected pullback metric is positive definite on
the full one-dimensional path. A robust Krawczyk continuation has contraction
ratio below `0.85`.

The proof has two distinct combinatorial layers. The metric layer has six
patches, five bridges, and 1003 gap-free interval segments. The robust root
layer has seven patches and seven bridges. Exact rational replay establishes
positivity of the three serialized leading-principal-minor bounds on every
metric segment. A separate implementation reconstructs the formula on the
weakest-third-minor segment selected from each metric patch and obtains 18
strictly positive interval minors.

At the endpoint `d=1.5`, an independently constructed 83-real Krawczyk box is
nested with the continuation endpoint. Independent baseline and correction
metric implementations, followed by exact Fraction determinant bounds, prove
local positivity. A fixed rational direction ensemble proves existence of a
relative correction exceeding 1.83%.

## Evidence boundary

The six stratified formula reconstructions reuse certified rational root
boxes. The remaining 997 segments have exact replay of serialized bounds but
not a second formula implementation. The theorem concerns a one-dimensional
trajectory, not an open six-real-dimensional tube or the full compact
hypersurface. It supplies no global Ricci/Monge--Ampere error bound and no
complete physical-vacuum theorem.

## Reproducibility boundary

The publication traversal has zero unresolved and zero ambiguous SHA-256
edges. Three unsuccessful diagnostics are excluded from proof traversal and
retained byte-for-byte in a separate historical manifest. The repository's
stdlib-only checker verifies the release manifest and recomputes the exact
serialized inequalities; heavy Arb generation remains upstream.
