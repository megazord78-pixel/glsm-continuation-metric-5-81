# Certificate schema

All proof-bearing real numbers used by the portable verifier are serialized
as reduced rational pairs:

```json
{"numerator": 29997, "denominator": 20000}
```

The denominator must be positive and the pair must already be reduced.  The
checker rejects noncanonical pairs.

## Central certificate

`results/2026-08-13_string_5_81_glsm_continuation_chain_publication_certificate.json`
is the human-readable summary produced by the checker.  Its main fields are:

- `proof_data.certified_depth_interval`: exact connected interval;
- `proof_data.minimum_leading_principal_minor_lower_bounds`: global minima
  over the 1003 metric segments;
- `proof_data.robust_root_chain`: seven-patch/seven-bridge Krawczyk replay;
- `proof_data.independent_local_material_box`: endpoint root, determinants,
  and correction-materiality witness;
- `proof_data.stratified_independent_metric_formula_replication`: six
  independently reconstructed weakest segments;
- `proof_data.transitive_provenance_audit`: resolved proof graph;
- `status`: admitted and explicitly rejected claim classes.

## Metric covers

Each segment contains exact lower and upper parameter endpoints and three
leading-principal-minor records.  The verifier requires dimensions
`[1, 2, 3]`, strict positive lower bounds, consecutive indices, exact
adjacency within a patch, and connected union equal to
`[29997/20000, 3/2]`.

## Krawczyk records

Every bridge or patch carries 83 componentwise image upper bounds and target
radii.  The only admitted relation is strict inclusion `image < radius`.
The robust chain additionally requires every ratio to be strictly below
`17/20`.

## Endpoint metric

The baseline and correction are `3 x 3` complex interval matrices.  Each
entry has directed rational real and imaginary bounds.  The checker combines
the matrices and recomputes the first three determinants by the Leibniz
formula using exact rational interval arithmetic.  Serialized determinant
bounds are accepted only when they agree exactly with this recomputation.

## Manifest

`MANIFEST.sha256.json` lists every release file except itself.  Generated
caches, bytecode, and build products are not part of the release.

