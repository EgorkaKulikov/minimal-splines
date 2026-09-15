package splines

import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import splines.ReferenceSplines
import splines.functionals.AveragingFunctionals
import splines.functionals.DeBoorFixFunctionals
import splines.functionals.FunctionalFamily
import splines.functionals.ProjFunctionals
import splines.functionals.ThreePointFunctionals
import splines.nonDegenerate
import org.junit.jupiter.api.Tag
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertTrue

/**
 * Invariants of the basis and the functionals: partition of unity, biorthogonality, exactness on
 * span φ.
 *
 * Every check runs on four grids — uniform, quasi-uniform, graded and geometric — in order to
 * catch errors that show up only for unequal steps.
 */
@Tag("fast")
class SplineCoreHealthCheckTest {

    private companion object {
        /** Comparison threshold for identities that hold exactly (up to rounding). */
        const val EXACT_IDENTITY_TOLERANCE = 1e-10

        /** Threshold for quantities that accumulate the error of solving small linear systems. */
        const val LINEAR_SOLVE_TOLERANCE = 1e-9

        /** Threshold for quantities that accumulate the error of projection and spline evaluation. */
        const val PROJECTION_TOLERANCE = 1e-8

        /**
         * Lower bound for the defect of the closed formula on a non-uniform grid: the value is
         * deliberately below the actual defect (about 0.2), so the check
         * [closedFormIsUniformGridOnly] does not depend on the exact value.
         */
        const val CLOSED_FORM_GRADED_DEFECT_FLOOR = 1e-2

        /** Number of sample points inside the interval used in pointwise comparisons. */
        const val SAMPLE_COUNT = 200
    }

    private val uniformGrid = Grid.uniform(8)
    private val quasiUniformGrid = Grid.quasiUniform(8)
    private val gradedGrid = Grid.graded(8)
    private val geometricGrid = Grid.geometric(8)
    private val sampleFractions = (0..SAMPLE_COUNT).map { it.toDouble() / SAMPLE_COUNT }
    private val allSystems = listOf(GeneratingSystem.B, GeneratingSystem.H, GeneratingSystem.T)

    /**
     * The test grids together with names used in failure messages.
     *
     * The `graded` and `geometric` grids are substantially non-uniform: for `graded` the ratio of
     * neighbouring steps is fixed for any n, and it is exactly this grid that is used when computing
     * the non-uniform tables (`verification.Sec4VerificationTool`).
     */
    private val namedGrids: List<Pair<String, Grid>> = listOf(
        "uniform" to uniformGrid,
        "quasiUniform" to quasiUniformGrid,
        "graded" to gradedGrid,
        "geometric" to geometricGrid,
    )

    /** Returns the largest deviation over all test grids. */
    private fun worstOverGrids(action: (Grid) -> Double): Double =
        namedGrids.maxOf { (_, grid) -> action(grid) }

    /** Sample points uniformly covering the grid interval. */
    private fun samplePoints(grid: Grid): List<Double> =
        sampleFractions.map { grid.a + (grid.b - grid.a) * it }

    /**
     * The general minimal-spline basis on the polynomial generating system must coincide with the
     * classical explicit formula for the quadratic B-spline — both by value and by the first
     * derivative.
     *
     * This is the key check: the general construction via inversion of the matrix `M_k` is verified
     * against an independently written closed formula.
     */
    @Test
    fun polynomialSplineMatchesClosedFormB() {
        val deviation = worstOverGrids { grid ->
            val basis = MinimalSplineBasis(GeneratingSystem.B, grid)
            var worst = 0.0
            for (t in samplePoints(grid)) {
                for (j in -2..grid.n - 1) if (nonDegenerate(grid, j)) {
                    worst = maxOf(worst, abs(basis.omega(j, t) - ReferenceSplines.omegaB(grid, j, t)))
                    worst = maxOf(
                        worst,
                        abs(basis.omegaDeriv(j, t) - ReferenceSplines.omegaBDeriv(grid, j, t)),
                    )
                }
            }
            worst
        }
        assertTrue(
            deviation < EXACT_IDENTITY_TOLERANCE,
            "The B basis must coincide with the explicit quadratic B-spline formula, " +
                "largest deviation = $deviation",
        )
    }

    /**
     * The same for the hyperbolic generating system: the general basis must coincide with the
     * explicit formula for the hyperbolic minimal spline.
     */
    @Test
    fun hyperbolicSplineMatchesClosedFormH() {
        val deviation = worstOverGrids { grid ->
            val basis = MinimalSplineBasis(GeneratingSystem.H, grid)
            var worst = 0.0
            for (t in samplePoints(grid)) {
                for (j in -2..grid.n - 1) if (nonDegenerate(grid, j)) {
                    worst = maxOf(worst, abs(basis.omega(j, t) - ReferenceSplines.omegaH(grid, j, t)))
                }
            }
            worst
        }
        assertTrue(
            deviation < EXACT_IDENTITY_TOLERANCE,
            "The H basis must coincide with the explicit hyperbolic spline formula, " +
                "largest deviation = $deviation",
        )
    }

    /**
     * Partition of unity: the sum of all basis splines is identically equal to one at any point of
     * the interval. The property follows from the first component of the generating vector function
     * being a constant.
     */
    @Test
    fun basisFormsPartitionOfUnity() {
        val deviation = worstOverGrids { grid ->
            var worst = 0.0
            for (system in allSystems) {
                val basis = MinimalSplineBasis(system, grid)
                for (t in samplePoints(grid)) {
                    var sum = 0.0
                    for (j in -2..grid.n - 1) sum += basis.omega(j, t)
                    worst = maxOf(worst, abs(sum - 1.0))
                }
            }
            worst
        }
        assertTrue(
            deviation < EXACT_IDENTITY_TOLERANCE,
            "The sum of the basis splines must equal one, largest deviation = $deviation",
        )
    }

    /**
     * Biorthogonality of the projection functionals: `theta_i(omega_j) = delta_ij`.
     * It is precisely this property that makes the operator `P_theta` a projector.
     */
    @Test
    fun projectionFunctionalsAreBiorthogonal() {
        val deviation = worstOverGrids { grid ->
            var worst = 0.0
            for (system in allSystems) {
                val basis = MinimalSplineBasis(system, grid)
                val funcs = ProjFunctionals(basis)
                for (i in -2..grid.n - 1) for (j in -2..grid.n - 1) {
                    val value = funcs.chi(i).apply(
                        { t -> basis.omega(j, t) },
                        { t -> basis.omegaDeriv(j, t) },
                    )
                    worst = maxOf(worst, abs(value - if (i == j) 1.0 else 0.0))
                }
            }
            worst
        }
        assertTrue(
            deviation < LINEAR_SOLVE_TOLERANCE,
            "The theta functionals must be biorthogonal to the basis, largest deviation = $deviation",
        )
    }

    /**
     * Biorthogonality of the de Boor–Fix functionals for all three variants `r = 0, 1, 2`,
     * including the boundary indices.
     *
     * For the boundary functionals (`j = -2` and `j = n-1`) the implementation uses plain values at
     * the endpoints of the interval — an assumption not written out explicitly in the primary
     * source, and the present check serves as its justification.
     */
    @Test
    fun deBoorFixFunctionalsAreBiorthogonalForAllOrders() {
        val deviation = worstOverGrids { grid ->
            var worst = 0.0
            for (system in allSystems) {
                val basis = MinimalSplineBasis(system, grid)
                for (r in 0..2) {
                    val funcs = DeBoorFixFunctionals(basis, r)
                    for (i in -2..grid.n - 1) for (j in -2..grid.n - 1) {
                        val value = funcs.chi(i).apply(
                            { t -> basis.omega(j, t) },
                            { t -> basis.omegaDeriv(j, t) },
                            { t -> basis.omegaDeriv2(j, t) },
                        )
                        worst = maxOf(worst, abs(value - if (i == j) 1.0 else 0.0))
                    }
                }
            }
            worst
        }
        assertTrue(
            deviation < LINEAR_SOLVE_TOLERANCE,
            "The functionals xi<0>, xi<1>, xi<2> must be biorthogonal to the basis, " +
                "largest deviation = $deviation",
        )
    }

    /**
     * Idempotence of the projectors: if the function is already a spline, its projection must return
     * exactly the same coefficients (`P^2 = P`).
     *
     * Checked only for the projectors (theta and all variants of xi). The families `mu` and `lambda`
     * are quasi-interpolants and are not required to be idempotent.
     *
     * The coefficients are pseudo-random with a FIXED seed: the test must be reproducible.
     */
    @Test
    fun projectorsAreIdempotentOnSplines() {
        val random = kotlin.random.Random(seed = 777)
        val deviation = worstOverGrids { grid ->
            var worst = 0.0
            for (system in allSystems) {
                val basis = MinimalSplineBasis(system, grid)
                val projectorFamilies = listOf(
                    ProjFunctionals(basis),
                    DeBoorFixFunctionals(basis, 0),
                    DeBoorFixFunctionals(basis, 1),
                    DeBoorFixFunctionals(basis, 2),
                )
                for (funcs in projectorFamilies) {
                    val coeffs = DoubleArray(grid.n + 2) { random.nextDouble(-1.0, 1.0) }
                    val spline = { t: Double -> basis.evalSpline(coeffs, t) }
                    val splineDeriv = { t: Double -> basis.evalSplineDeriv(coeffs, t) }
                    val splineDeriv2 = { t: Double -> basis.evalSplineDeriv2(coeffs, t) }
                    val projected = funcs.projectorCoeffs(spline, splineDeriv, splineDeriv2)
                    for (i in coeffs.indices) worst = maxOf(worst, abs(projected[i] - coeffs[i]))
                }
            }
            worst
        }
        assertTrue(
            deviation < PROJECTION_TOLERANCE,
            "The projectors theta and xi must be idempotent on splines, " +
                "largest deviation = $deviation",
        )
    }

    /**
     * Exactness on the generating space: each of the four functional families must reproduce
     * functions from `span{1, rho, sigma}` with no error.
     *
     * This is the minimal requirement on an approximation operator; its violation means an error in
     * the construction of the functional coefficients.
     *
     * The whole span is checked: all three generators `1`, `rho`, `sigma` and their linear
     * combination. The latter detects errors that cancel on the individual generators but not on a
     * general element of the space.
     *
     * The coefficients of the combination are pseudo-random with a fixed seed — the test must be
     * reproducible (the same convention as in [projectorsAreIdempotentOnSplines]).
     *
     * Derivatives up to the second order are passed: the xi family with `r = 0` uses `fDD`, and
     * without the third argument it would silently receive a zero second derivative. The family list
     * uses `r = 1`; adding `r = 0` must not weaken the check.
     */
    @Test
    fun allFamiliesAreExactOnGeneratingSpan() {
        val random = kotlin.random.Random(seed = 20240117)
        var deviation = 0.0
        var worstLabel = ""
        for ((gridName, grid) in namedGrids) {
            for (system in allSystems) {
                val basis = MinimalSplineBasis(system, grid)
                val families: List<FunctionalFamily> = listOf(
                    ProjFunctionals(basis),
                    DeBoorFixFunctionals(basis),
                    AveragingFunctionals(basis),
                    ThreePointFunctionals(basis),
                )
                val c0 = random.nextDouble(-1.0, 1.0)
                val c1 = random.nextDouble(-1.0, 1.0)
                val c2 = random.nextDouble(-1.0, 1.0)
                val members = listOf(
                    SpanMember("1", { 1.0 }, { 0.0 }, { 0.0 }),
                    SpanMember("rho", system.rho, system.rhoD, system.rhoDD),
                    SpanMember("sigma", system.sigma, system.sigmaD, system.sigmaDD),
                    SpanMember(
                        "combination",
                        { t -> c0 + c1 * system.rho(t) + c2 * system.sigma(t) },
                        { t -> c1 * system.rhoD(t) + c2 * system.sigmaD(t) },
                        { t -> c1 * system.rhoDD(t) + c2 * system.sigmaDD(t) },
                    ),
                )
                for (funcs in families) {
                    for (member in members) {
                        val projected = funcs.projectorCoeffs(member.g, member.gD, member.gDD)
                        for (t in samplePoints(grid)) {
                            val error = abs(member.g(t) - basis.evalSpline(projected, t))
                            if (error > deviation) {
                                deviation = error
                                worstLabel = "${funcs.name}/${system.name}/$gridName/${member.label} at t=$t"
                            }
                        }
                    }
                }
            }
        }
        assertTrue(
            deviation < PROJECTION_TOLERANCE,
            "All functional families must be exact on span{1, rho, sigma}, " +
                "largest deviation = $deviation ($worstLabel)",
        )
    }

    /**
     * An element of the generating space together with its two derivatives and a name for the
     * failure message. The derivatives are kept next to the function itself because the xi family
     * requires them simultaneously with the value, and a mismatch between them would produce not a
     * compilation error but a silently wrong check.
     */
    private data class SpanMember(
        val label: String,
        val g: (Double) -> Double,
        val gD: (Double) -> Double,
        val gDD: (Double) -> Double,
    )

    /**
     * Check of the closed formula for the projection functionals: on a uniform grid with the
     * polynomial generating system the coefficients of the interior functionals must equal
     * `{1/14, -2/7, 10/7, -2/7, 1/14}`.
     *
     * The numbers are taken from the primary source (see `docs/REFERENCES.md`, section
     * “Approximation functionals”), so the check verifies the implementation against the published
     * formula rather than against itself.
     */
    @Test
    fun closedFormCoefficientsMatchPublishedValues() {
        val grid = Grid.uniform(8)
        val basis = MinimalSplineBasis(GeneratingSystem.B, grid)
        val funcs = ProjFunctionals(basis)
        val expected = doubleArrayOf(1.0 / 14.0, -2.0 / 7.0, 10.0 / 7.0, -2.0 / 7.0, 1.0 / 14.0)
        var deviation = 0.0
        for (j in 0..grid.n - 3) {
            val actual = funcs.closedFormInternal(j).coeffs
            for (k in actual.indices) deviation = maxOf(deviation, abs(actual[k] - expected[k]))
        }
        assertTrue(
            deviation < EXACT_IDENTITY_TOLERANCE,
            "The coefficients of the closed formula for theta must match the published " +
                "{1/14, -2/7, 10/7, -2/7, 1/14}, largest deviation = $deviation",
        )
    }

    /**
     * The closed formula is applicable ONLY to a uniform grid.
     *
     * The symmetric set of weights `{E^2, -DE, CD-BE, -DE, E^2}/K1` yields biorthogonality only when
     * `A = E` and `B = D` — this is the symmetry of the coordinate spline about the centre of its
     * support, which holds on a uniform grid for an even/odd pair of generators. Outside this case
     * the formula does not reproduce even a constant, so it is used only as a cross-check on a
     * uniform grid (see [closedFormCoefficientsMatchPublishedValues]), while the working path
     * remains the solution of the local biorthogonality system.
     *
     * Reproduction of a constant is checked: the value of the functional on `f = 1` equals the sum
     * of its weights, since `sum_j omega_j = 1`.
     */
    @Test
    fun closedFormIsUniformGridOnly() {
        val one = { _: Double -> 1.0 }
        val zero = { _: Double -> 0.0 }
        var builtDefect = 0.0
        for (system in allSystems) {
            val basis = MinimalSplineBasis(system, gradedGrid)
            val funcs = ProjFunctionals(basis)
            var closedDefect = 0.0
            for (j in 0..gradedGrid.n - 3) {
                builtDefect = maxOf(builtDefect, abs(funcs.chi(j).apply(one, zero) - 1.0))
                closedDefect = maxOf(closedDefect, abs(funcs.closedFormInternal(j).apply(one, zero) - 1.0))
            }
            // Actual defects at n = 8, ratio = 2: B — 0.2222222222222296,
            // H — 0.2208866259181210, T — 0.2235658323679606. The exact value for the
            // polynomial generating system is 2/9 and does not decrease under refinement.
            assertTrue(
                closedDefect > CLOSED_FORM_GRADED_DEFECT_FLOOR,
                "The closed formula is not biorthogonal on a non-uniform grid, hence the defect of " +
                    "reproducing a constant must be large; system ${system.name}, " +
                    "defect = $closedDefect",
            )
        }
        // Actual defect of the main path at n = 8, ratio = 2: 1.7e-13.
        assertTrue(
            builtDefect < LINEAR_SOLVE_TOLERANCE,
            "The main construction path for theta must reproduce a constant on any grid, " +
                "largest deviation = $builtDefect",
        )
    }
}
