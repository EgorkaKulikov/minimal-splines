package splines.metrics

import splines.Grid
import org.junit.jupiter.api.Tag
import kotlin.math.abs
import kotlin.math.exp
import kotlin.math.sin
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Tests of the control-grid parameter of [errorEh]: measuring the uniform norm on a refined grid.
 */
@Tag("fast")
class MetricsControlGridTest {

    /** Refining the control grid can only increase the maximum over the points, never decrease it. */
    @Test fun finerControlGridNeverDecreasesMaximum() {
        val g = Grid.uniform(16, 0.0, 1.0)
        val exact = { t: Double -> t * t }
        val eval = { t: Double -> t * t + 0.01 * sin(50.0 * t) }
        val coarse = errorEh(exact, eval, g)
        val fine = errorEh(exact, eval, g, refinement = 1000)
        assertTrue(fine >= coarse - 1e-15, "coarse=$coarse, fine=$fine")
    }

    /**
     * For a smooth difference, going from 100n+1 to 1000n+1 points shifts the value by at most
     * 4.071e-4 in relative terms, i.e. the default control grid does not underestimate the error
     * maximum.
     */
    @Test fun defaultControlGridDoesNotUnderestimateSmoothError() {
        val g = Grid.uniform(32, 0.0, 1.0)
        val exact = { t: Double -> exp(-t) }
        val eval = { t: Double -> exp(-t) + 1e-4 * sin(9.0 * t) * (1.0 - t) }
        val coarse = errorEh(exact, eval, g)
        val fine = errorEh(exact, eval, g, refinement = 1000)
        val shift = abs(fine - coarse) / fine
        assertTrue(shift <= 4.071e-4, "relative shift $shift exceeded the reference 4.071e-4")
    }

    /**
     * Why the parameter is needed at all: for a difference with a narrow spike between the control
     * grid points, a coarse grid UNDERESTIMATES the maximum, and without refinement this stays unnoticed.
     */
    @Test fun coarseControlGridMissesNarrowSpike() {
        val g = Grid.uniform(4, 0.0, 1.0)
        val spikeCentre = 0.5 + 1.0 / (2.0 * 100 * 4) // exactly between the points of the 100n+1 grid
        val exact = { _: Double -> 0.0 }
        val eval = { t: Double -> exp(-4e7 * (t - spikeCentre) * (t - spikeCentre)) }
        val coarse = errorEh(exact, eval, g)
        val fine = errorEh(exact, eval, g, refinement = 100_000)
        assertTrue(coarse < 0.9, "the coarse control grid must miss the spike, got $coarse")
        assertTrue(fine > 0.99, "the refined control grid must detect the spike, got $fine")
    }

    /** refinement < 1 would give m = 0 and a 0/0 division: that is a contract violation, not a silent NaN. */
    @Test fun refinementBelowOneIsRejected() {
        val g = Grid.uniform(4)
        assertFailsWith<IllegalArgumentException> { errorEh({ 0.0 }, { 0.0 }, g, refinement = 0) }
        assertFailsWith<IllegalArgumentException> { errorEh({ 0.0 }, { 0.0 }, g, refinement = -1) }
    }

    /** refinement = 1 is legal: the control grid coincides with the nodes of the uniform base grid. */
    @Test fun refinementOneEvaluatesAtGridNodes() {
        val g = Grid.uniform(4, 0.0, 1.0)
        val e = errorEh({ _ -> 0.0 }, { t -> t }, g, refinement = 1)
        assertEquals(1.0, e, 1e-12)
    }
}
