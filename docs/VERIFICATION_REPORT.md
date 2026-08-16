# Verification report for release candidate 2.0.0

Date: 2026-08-16
Status: scientific verification gate passed.

## Completed checks

The dependency-free release gate passed in the source package and in an
independent copy located outside the publication tree. Both runs authenticated
all 72 manifested files and reported:

- 1,084 exact parent root boxes;
- 34,688 metric child boxes;
- 89,972 strict root-component inclusions;
- 1,083 positive longitudinal overlaps;
- exact 32-child transverse tiling at every path center;
- positive lower bounds for all 104,064 metric-minor intervals;
- independent exact replay of the endpoint root, metric, and materiality data;
- 11/11 positive and fail-closed tests.

The exact global metric-minor lower bounds were

```text
3459225413190877/288230376151711744
6861360338265545/18446744073709551616
2100827914296867/18889465931478580854784
```

The minimum root inclusion margin was

```text
9335426793689657913651/10910555522022028301723238400000
```

The exact-integer registry was reconstructed from the raw adapter with a
fail-closed decoder. A deliberately nonintegral binary64 mutation was rejected.

## Formula-level reconstruction

Using the pinned python-flint 0.9.0 environment, the preliminary smoke check
reconstructed metric child sequence 0 from the standalone raw inputs. Every
proof-relevant field was exactly equal to the archived record, including the
third-minor lower endpoint

```text
8461319195960853/75557863725914323419136
```

The diagnostic root ratio was also equal, although diagnostic floats are not
part of the proof comparison.

## Manuscript QA

The combined path-plus-tube manuscript compiled to an eight-page A4 PDF after
BibTeX and two final
pdfLaTeX passes. The final compiler pass contained no LaTeX warnings or
overfull boxes. All eight pages were rendered to PNG and visually inspected;
no clipped text, overlaps, broken glyphs, or unreadable table entries were
observed.

## Full regeneration result

A fresh eight-worker run regenerated all 34,688 metric children and all 1,084
parent roots under python-flint 0.9.0. Every proof-relevant regenerated field
was exactly equal to the corresponding field in both archived compressed
certificates. Diagnostic floating-point fields were deliberately excluded from
the equality relation and were not used as proof inputs.

The full release gate reported:

```text
full_arb_regeneration_passed=true
proof_fields_equal_to_release=true
metric_child_boxes=34688
parent_root_boxes=1084
```

Scientific verification no longer blocks the publication workflow. Git
tagging, archival release creation, and DOI metadata remain publication
operations rather than proof obligations.
