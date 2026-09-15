# Golden references of the library behaviour

The `*.json` files record the numerical behaviour of `minimal-splines` version 0.1.0 and serve
as regression control during rework: `grid.json` (grid nodes and flags),
`basis.json` (basis and spline at the control points), `functionals.json` (families of
approximation functionals), `metrics.json` (`errorEh`),
`algebra.json` (`nonDegenerate`, generating systems).

Format: JSON without external dependencies; every `Double` is written as 16 hexadecimal
digits of raw bits (`Double.toRawBits`), which rules out losses during printing and parsing. The
`generatedWith` field holds the library version, the name of the BLAS/LAPACK implementation and the generation date.

Regeneration: `./gradlew regenerateGolden` (the task is excluded from `test`, `fastTest` and Kover).

Verification tolerances (`Golden*Test`, tag `fast`): a bit-for-bit match for integer and
combinatorial logic (`interval`, `isCoincident`, `nonDegenerate`, family flags) and for explicit
formulas (`Grid`, algebra of 3-vectors, `phi`, `wronskian` of systems without transcendental functions); the relative
`max|got - exp| / max(max|exp|, 1) <= 1e-10` for quantities that go through the solution of small linear systems
(`omega*`, `evalSpline*`, projector coefficients, `cChi`, `errorEh`), where a discrepancy in the
order of operations of the BLAS/LAPACK implementation and in the way 3x3 matrices are inverted is admissible.

Decision 1c: for `basis`, `functionals` and `errorEh` the relative tolerance equals `1e-10`
(`GoldenCompare.BASIS_TOLERANCE`). The 0.1.0 golden references were produced by inverting M_k with Cramer's rule in the
global coordinates of the generating system, at condition numbers up to ~10⁴; the error of that
path is of order cond·ε ≈ 10⁻¹²…10⁻¹¹, which is exactly what the comparison with LAPACK shows (max. 1.05·10⁻¹¹).
A tolerance of 10⁻¹⁰ covers the golden reference error with a tenfold margin and is orders of magnitude smaller than any
method error. For `grid` the comparison is bit for bit.

Quantities computed through `Math.sin`, `Math.cos`, `Math.pow` (generating system T, grids
`quasiUniform` and `geometric`) are compared with the relative tolerance
`GoldenCompare.TRANSCENDENTAL_TOLERANCE = 16 ulp(1) ≈ 3.6e-15` rather than bit for bit: HotSpot uses
platform intrinsics for these functions, whose results differ in the last bit between x86_64 and aarch64
(the golden references were taken on macOS aarch64, while CI also runs on ubuntu x86_64). In `geometric` a one-ulp
discrepancy in q = R^{1/(n-1)} is carried over to q^j with a factor of j ≤ n and adds up with the discrepancies
of the individual powers — up to ~7 ulp(1) at n = 32; a tolerance of 16 ulp(1) covers this with a margin and stays
orders of magnitude below any algorithmic discrepancy.
