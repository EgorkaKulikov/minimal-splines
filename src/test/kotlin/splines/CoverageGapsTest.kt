package splines

import org.junit.jupiter.api.Tag
import splines.functionals.DeBoorFixFunctionals
import splines.functionals.DiscreteDeBoorFixFunctionals
import splines.metrics.errorEh
import kotlin.math.sinh
import kotlin.test.Test
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Rejection branches for invalid input and numerical degeneracy that are not reached by the main
 * tests: the parameters of the de Boor--Fix families, the parameters of the grid factories and of
 * the metrics, degeneracy of the approximation relation for a user-defined system in global
 * coordinates, overflow of the global representation of H and degeneracy of xi^<0> at an inflection
 * point.
 */
@Tag("fast")
class CoverageGapsTest {
    private val basis = MinimalSplineBasis(GeneratingSystem.B, Grid.uniform(8, 0.0, 1.0))

    private fun message(e: Throwable): String = e.message ?: ""

    @Test
    fun deBoorFix_orderOutsideRangeRejected() {
        for (r in listOf(-1, 3)) {
            val e = assertFailsWith<IllegalArgumentException>("DeBoorFix r=$r") { DeBoorFixFunctionals(basis, r = r) }
            assertTrue(message(e).contains("parameter r must be 0, 1 or 2") && message(e).contains("got $r"), message(e))
        }
        for (r in listOf(0, 3)) {
            val e = assertFailsWith<IllegalArgumentException>("DiscreteDeBoorFix r=$r") { DiscreteDeBoorFixFunctionals(basis, r = r) }
            assertTrue(message(e).contains("parameter r must be 1 or 2") && message(e).contains("got $r"), message(e))
        }
    }

    @Test
    fun gridFactories_invalidRatioRejected() {
        val geometric = assertFailsWith<IllegalArgumentException> { Grid.geometric(8, R = 0.0) }
        assertTrue(message(geometric).contains("R must be > 0"), message(geometric))
        val negative = assertFailsWith<IllegalArgumentException> { Grid.geometric(8, R = -2.0) }
        assertTrue(message(negative).contains("R=-2.0"), message(negative))
    }

    @Test
    fun errorEh_invalidRefinementRejected() {
        val e = assertFailsWith<IllegalArgumentException> { errorEh({ 0.0 }, { 0.0 }, basis.grid, refinement = 0) }
        assertTrue(message(e).contains("refinement must be at least 1"), message(e))
    }

    /**
     * A user-defined system without a local representation is evaluated in global coordinates;
     * for (1, t, t^2) on [0, 1e6] the condition number of T_k M_k exceeds MAX_CONDITION, and the
     * constructor rejects the grid with a message about the conditioning of the approximation
     * relation matrix.
     */
    @Test
    fun customSystem_illConditionedGlobalMatrixRejected() {
        val globalPolynomial = GeneratingSystem(
            "poly-global", rho = { t -> t }, sigma = { t -> t * t },
            rhoD = { 1.0 }, sigmaD = { t -> 2.0 * t },
            rhoDD = { 0.0 }, sigmaDD = { 2.0 },
        )
        val e = assertFailsWith<IllegalArgumentException> { MinimalSplineBasis(globalPolynomial, Grid.uniform(100, 0.0, 1e6)) }
        assertTrue(
            message(e).contains("Approximation relation matrix on interval") && message(e).contains("ill-conditioned"),
            message(e),
        )
    }

    /**
     * The basis H on [0, 2000] with h = 1 is built in local coordinates, but the global representation
     * sinh(t), cosh(t) overflows for t > 710: the diagnostic computeA(j) rejects such a j.
     */
    @Test
    fun globalComputeA_overflowRejected() {
        val hyperbolic = MinimalSplineBasis(GeneratingSystem.H, Grid.uniform(2000, 0.0, 2000.0))
        assertTrue(sinh(1501.0).isInfinite())
        val e = assertFailsWith<IllegalArgumentException> { hyperbolic.computeA(1500) }
        assertTrue(message(e).contains("computeA(j=1500)") && message(e).contains("overflow (scale"), message(e))
        // Near the start of the interval the global representation is finite and has the same size as the local one.
        assertTrue(hyperbolic.computeA(0).size == 3)
    }

    /**
     * The functionals xi^<0> require rho' sigma'' - rho'' sigma' != 0; for the system (1, t, t^3) this
     * quantity vanishes at the inflection point t = 0, and with a grid node at zero the family is not
     * built.
     */
    @Test
    fun deBoorFixZero_inflectionNodeRejected() {
        val cubic = GeneratingSystem(
            "cubic", rho = { t -> t }, sigma = { t -> t * t * t },
            rhoD = { 1.0 }, sigmaD = { t -> 3.0 * t * t },
            rhoDD = { 0.0 }, sigmaDD = { t -> 6.0 * t },
        )
        val cubicBasis = MinimalSplineBasis(cubic, Grid.uniform(8, -1.0, 1.0))
        val e = assertFailsWith<IllegalArgumentException> { DeBoorFixFunctionals(cubicBasis, r = 0) }
        assertTrue(message(e).contains("degenerate denominator rho'sigma''-rho''sigma'="), message(e))
    }
}
