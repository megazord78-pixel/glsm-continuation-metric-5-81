# Publication checklist

## Before public release

- [x] Freeze the theorem and evidence boundary.
- [x] Package the exact interval records and a dependency-free verifier.
- [x] Include negative mutation tests.
- [x] Include a complete LaTeX manuscript and bibliography.
- [x] Include BSD-3-Clause licensing and `CITATION.cff`.
- [x] Bind all files with a release manifest.
- [x] Compile `manuscript/main.tex` with BibTeX and visually inspect all PDF pages.
- [x] Run `python verify_release.py` on a clean unpack of the final ZIP.
- [x] Create a dedicated public GitHub repository from this directory only.
- [x] Add GitHub Actions for Python 3.11, 3.12, and 3.13.
- [x] Tag the initial immutable commit as `v1.0.0`.
- [x] Archive the metadata-only `v1.0.1` update with Zenodo and insert DOI
  `10.5281/zenodo.21923926` into `CITATION.cff`.

## Preprint submission

- [ ] Upload `main.tex`, `references.bib`, and the generated `.bbl` to arXiv.
- [ ] Use the title and abstract from `manuscript/main.tex` unchanged.
- [ ] Select primary category `hep-th`; consider `math.NA` as a cross-list
  only if the endorsement/category rules allow it.
- [ ] Put the GitHub release and Zenodo DOI in the comments or data-availability
  field.
- [ ] Check that no sentence promotes the theorem to a global CY metric or a
  physical alpha-prime correction.

## Journal package

- [ ] Prepare a short cover letter emphasizing validated numerics, the
  explicit `(5,81)` geometry, and the narrow theorem boundary.
- [ ] Suggest reviewers with expertise in numerical Calabi--Yau metrics and
  validated numerics; avoid anyone with a conflict of interest.
- [ ] Submit to a no-mandatory-APC route.  Suitable scope candidates should be
  checked at submission time because policies can change.
- [ ] State that code and exact proof records are permanently archived.

## Priority language

Use only:

> To the best of our knowledge, this is the first computer-assisted
> interval-arithmetic positivity theorem for a nontrivially corrected GLSM
> quotient metric along a connected continuation path on an explicit compact
> Calabi--Yau threefold used in a string-compactification model.

Do not shorten this to “the first corrected Calabi--Yau metric theorem.”
