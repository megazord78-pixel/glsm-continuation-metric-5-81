# Certificate schema

`data/root_certificate.json.gz` contains 1,084 parent records. Each record
contains an exact rational centre/radius box, directed hypersurface Newton
bounds, 81 directed D-term Krawczyk bounds, and diagnostic-only float history.

`data/metric_certificate.json.gz` contains 34,688 child records. Each contains
three directed rational intervals for the leading principal minors. The 32
children per parent exactly tile the transverse cube.

`data/certificate.json` stores only summaries independently recomputed by the
quick verifier. Its Boolean scope fields restrict rather than establish the
theorem.

Every rational pair is reduced and encoded as

```json
{"numerator": 1, "denominator": 100000000}
```

Every directed interval is encoded as

```json
{"lower": {"numerator": 1, "denominator": 2},
 "upper": {"numerator": 2, "denominator": 3}}
```
