# Trust model

## Certified statement

The release certifies one connected branch over

\[
d\in[1.49985,1.5]
\]

and positive definiteness of the frozen corrected `3 x 3` Hermitian
pullback metric on that branch.  It also certifies an independently rebuilt
83-real local root and corrected metric at `d=1.5`, where at least one fixed
rational direction has a relative correction larger than one percent.

## Proof-bearing layers

1. Arb-generated interval records enclose the implicit hypersurface and
   81-component moment-map branch, its derivatives, and the metric entries.
2. Seven root patches and seven root bridges prove existence and local
   uniqueness of one connected branch with componentwise Krawczyk ratios
   below `0.85`.
3. Six metric patches contain 1003 exactly adjacent rational parameter
   segments.  Positive lower bounds for all three leading principal minors
   prove Hermitian positive definiteness on every segment.
4. A formula-independent implementation rebuilds the metric on the weakest
   third-minor segment of each patch.  It checks 18 positive interval minors,
   while reusing the certified root boxes.
5. At the endpoint, a separate 83-real Krawczyk construction and independent
   baseline/correction implementations recompute the three determinants with
   exact `fractions.Fraction` interval arithmetic.
6. The portable checker authenticates every packaged attachment and replays
   the final rational inequalities without NumPy, Arb, JAX, CYTools, or the
   parent project.

## Trusted base

The mathematical result depends on:

- the disclosed toric chart, charge matrix, hypersurface, FI data, and frozen
  correction ansatz used by the upstream generators;
- correctness of Arb/python-flint and the generator implementations that
  produced the directed rational enclosures;
- the serialized exact-rational proof records in this release;
- the complete frozen correction parameter vector in
  `data/frozen_correction_parameters.json`;
- Python integer and `fractions.Fraction` arithmetic in the portable checker;
- SHA-256 for binding the release files and provenance graph.

Approximate centers, floating-point radius proposals, fitted coefficients,
and numerical inverse matrices are candidate witnesses.  They become
proof-bearing only after the interval inclusions and rational inequalities
are checked.

## Deliberate independence boundaries

The portable checker does not reevaluate the analytic GLSM formula with Arb.
Conversely, the six stratified formula reconstructions do not independently
prove the reused root boxes.  These two boundaries are stated explicitly so
that neither serialized replay nor partial formula replication is presented
as a complete second implementation.

## Claims excluded from the release

- positivity on an open six-real-dimensional tube;
- positivity on the full compact Calabi--Yau threefold or all 389 charts;
- a global Ricci-flatness or Monge--Ampere residual bound;
- derivation of the frozen ansatz as the physical alpha-prime correction;
- stabilization of a complete string vacuum, uplift, or phenomenology;
- global or model-independent uniqueness of the continuation branch.

The word *corrected* always refers to the exact frozen `ddbar` ansatz stated
in the manuscript, not to a completed string-effective-action calculation.

## Fail-closed behavior

- A missing, added, or hash-mismatched release file fails verification.
- A gap in the parameter cover fails verification.
- A zero or negative principal-minor lower bound fails verification.
- Equality at a Krawczyk-box boundary fails verification.
- A stale endpoint or provenance hash fails verification.
- Mutation tests exercise each of these central failure modes.
