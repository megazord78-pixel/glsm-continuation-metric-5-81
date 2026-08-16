# Trust model

## Quick verification

The dependency-free verifier trusts Python's integer arithmetic, `fractions`,
JSON/gzip parsing and SHA-256. It independently derives all final strict
rational inequalities and product-cover relations. Stored Boolean verdicts and
float diagnostics are not proof inputs.

It does not derive the Arb enclosures from the GLSM formulas. Therefore it is a
last-mile certificate replay, not a full formula-level regeneration.

## Full regeneration

The optional full mode additionally trusts Python, python-flint/Arb/FLINT, the
frozen exact geometry registry, the disclosed formula implementation and the
operating environment. It rebuilds the parent root and metric child enclosures
and compares them with the release certificates.

## Common-mode boundary

Neither mode proves that the selected GLSM chart or corrected Kähler ansatz is
the physically complete effective theory. The model-to-equations appendix and
exact-integer registry make that boundary explicit; independent external
replication remains stronger evidence.
