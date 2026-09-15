package splines.metrics

import splines.Grid
import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Tests of the error metric errorEh (uniform norm on a fine control grid).
 * The convergence orders `orders` and the constant `constCh` live in `numerical-core`
 * (`numerics.ConvergenceRates`) and are tested there (`ConvergenceRatesTest`).
 */
@Tag("fast")
class MetricsTest {
    /** errorEh: for an exact/eval pair with a known difference it yields the exact maximum modulus. */
    @Test fun errorEhCapturesMaxDifference() {
        val g = Grid.uniform(4, 0.0, 1.0)
        // eval = exact + bump; max |bump| = 0.5 at t=1 (linear shift 0.5*t)
        val e = errorEh({ t -> t }, { t -> t + 0.5 * t }, g)
        assertEquals(0.5, e, 1e-12)
    }

    /** errorEh = 0 for identical functions. */
    @Test fun errorEhZeroWhenEqual() {
        val g = Grid.uniform(4)
        assertEquals(0.0, errorEh({ t -> t * t }, { t -> t * t }, g), 1e-15)
    }

    /**
     * errorEh explicitly requires n >= 1: at n = 0 the divisor m = 100n would become zero.
     * Constructing Grid(0, ...) is no longer possible, so both links of the contract are checked:
     * a grid with n = 0 is unreachable, and the metric itself is computed on the smallest valid n = 1.
     */
    @Test fun errorEhRequiresPositiveN() {
        assertFailsWith<IllegalArgumentException> { Grid(0, doubleArrayOf(0.0)) }
        val e = errorEh({ t -> t }, { t -> t + 0.25 }, Grid.uniform(1, 0.0, 1.0))
        assertEquals(0.25, e, 1e-12)
    }

}
