package splines.functionals

import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

/**
 * Tests of the boolean family flags: theta is a projector without a derivative, xi is a projector
 * with a derivative, mu/lambda are quasi-interpolants. They cover the property getters.
 */
@Tag("fast")
class FlagsTest {
    private val basis = MinimalSplineBasis(GeneratingSystem.B, Grid.uniform(8))

    /** theta (ProjFunctionals): isProjector=true, usesDerivative=false. */
    @Test fun thetaFlags() {
        val theta = ProjFunctionals(basis)
        assertTrue(theta.isProjector)
        assertFalse(theta.usesDerivative)
        assertEquals("theta", theta.name)
    }

    /** The family names match the Greek notation. */
    @Test fun familyNames() {
        assertEquals("xi", DeBoorFixFunctionals(basis).name)
        assertEquals("mu", AveragingFunctionals(basis).name)
        assertEquals("lambda", ThreePointFunctionals(basis).name)
    }
}
