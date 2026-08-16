# Model-to-equations map

The exact data-to-formula boundary used by the article is:

| Mathematical object | Frozen source | Reconstruction |
|---|---|---|
| GLSM charge matrix `Q` | `data/exact_integer_geometry.json` | exact integers, shape 81 x 85 |
| Hypersurface exponents and coefficients | same registry | ten exponent vectors and Arb complex coefficient text |
| Selected toric chart | exact chart indices plus `regeneration/data/atlas.json` | three free complex directions and one eliminated coordinate |
| FI target | `regeneration/data/corrected_kahler_param.dat` and nef decomposition | exact decimal strings mapped through the disclosed rational basis |
| Frozen correction | `data/correction_parameters.json` | analytic second-order interval jets |
| Implicit root | hypersurface equation plus 81 moment maps | 83-real parameterized Krawczyk operator |
| Corrected metric | quotient pullback plus `ddbar(phi)` | three leading Hermitian principal minors |

The theorem starts after these choices are frozen. It does not prove that the
ansatz is the canonical Ricci-flat metric or the complete physical effective
metric.
