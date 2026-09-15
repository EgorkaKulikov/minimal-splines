package splines.functionals

import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import splines.metrics.errorEh
import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Targeted tests of the edge branches: the ValueFunctional contract, both maxOf branches in
 * errorEh (error peak in the middle), the default values of the Grid factories.
 */
@Tag("fast")
class EdgeBranchTest {
    /** ValueFunctional requires nodes and coeffs to have the same length (the require=false branch). */
    @Test fun valueFunctionalRejectsSizeMismatch() {
        assertFailsWith<IllegalArgumentException> {
            ValueFunctional(doubleArrayOf(0.0, 1.0), doubleArrayOf(1.0))
        }
    }

    /** errorEh with the difference peaking in the middle of the interval exercises both maximum-update branches. */
    @Test fun errorEhPeakInMiddle() {
        val g = Grid.uniform(4, 0.0, 1.0)
        // |diff| is a "hat": it grows up to t=0.5 and then decreases -> maxOf both updates and does not update
        val e = errorEh({ _ -> 0.0 }, { t -> 0.5 - Math.abs(t - 0.5) }, g)
        assertEquals(0.5, e, 1e-12) // maximum at t=0.5
    }

    /** Grid.uniform with default values (a=0,b=1) — covers the synthetic $default. */
    @Test fun gridUniformDefaults() {
        val g = Grid.uniform(4)
        assertEquals(0.0, g.a, 1e-15)
        assertEquals(1.0, g.b, 1e-15)
    }

    /** Grid.quasiUniform with default values — covers the synthetic $default. */
    @Test fun gridQuasiUniformDefaults() {
        val g = Grid.quasiUniform(6)
        assertEquals(0.0, g.a, 1e-15)
        assertEquals(1.0, g.b, 1e-15)
        for (i in 0 until g.n) assertTrue(g.x(i + 1) > g.x(i))
    }

    /** projectorCoeffs with an explicit derivative gD exercises the optional-parameter branch. */
    @Test fun projectorCoeffsWithExplicitDerivative() {
        val basis = MinimalSplineBasis(GeneratingSystem.B, Grid.uniform(8))
        val xi = DeBoorFixFunctionals(basis)
        val c = xi.projectorCoeffs({ t -> t }, { _ -> 1.0 })
        // xi reproduces the linear f(t)=t exactly
        for (t in listOf(0.2, 0.6)) assertEquals(t, basis.evalSpline(c, t), 1e-8)
    }
}
