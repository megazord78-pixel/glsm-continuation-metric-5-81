# Reproducibility

## Quick exact replay

Requirements: CPython 3.11 or later; no third-party packages.

```text
python verify_release.py --quick
```

Expected central counts:

- 1,084 root parent boxes;
- 34,688 metric child boxes;
- 89,972 exact componentwise root inclusions;
- 104,064 positive metric-minor intervals;
- 1,083 longitudinal gluing edges;
- 9 fail-closed tests.

The quick replay should take seconds to tens of seconds and should not create
proof artifacts.

## Full Arb regeneration

Use CPython 3.12 and install the pinned dependency:

```text
python -m pip install -r requirements-full.txt
python verify_release.py --full
```

The reference scientific environment was python-flint 0.9.0, FLINT 3.6.0,
100 Arb decimal digits, and eight workers. Worker counts can be set with
`GLSM_FULL_CUBE_WORKERS` and `GLSM_EXACT_ROOT_WORKERS`. The full run is
resumable under `regeneration/work/`; that directory is scratch state and is
excluded from the release manifest.

The comparison ignores diagnostic floating-point ratios, elapsed time,
checkpoint paths, and provenance hashes that necessarily change when code is
relocated. It requires exact equality of all directed rational proof bounds,
centers, radii, indices, and pass/fail proof fields.

## Manuscript

The reference PDF was built with MiKTeX/pdfLaTeX and BibTeX:

```text
cd manuscript
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

Compiler log files are build products, not release inputs.
