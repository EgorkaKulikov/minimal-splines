package splines

import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Tests for the geometric (non-uniform) grid Grid.geometric (the sn-article family):
 * the endpoints a,b are exact and the nodes at the endpoints are triple; strict monotonicity;
 * h = max step; the ratio of the extreme steps h_{n-1}/h_0 ~ R; local quasi-uniformity
 * h_j/h_{j-1} in [mu^{-1}, mu] with mu=q (for this grid the ratio is constant and equal to q).
 */
@Tag("fast")
class GridGeometricTest {
    private val tol = 1e-12

    @Test fun endpointsExactAndTripleKnots() {
        val g = Grid.geometric(16, 0.0, 1.0, R = 2.0)
        assertEquals(0.0, g.a, tol)
        assertEquals(1.0, g.b, tol)
        assertEquals(0.0, g.x(0), tol)
        assertEquals(1.0, g.x(g.n), tol)
        assertEquals(g.a, g.x(-2), tol); assertEquals(g.a, g.x(-1), tol)
        assertEquals(g.b, g.x(g.n + 1), tol); assertEquals(g.b, g.x(g.n + 2), tol)
    }

    @Test fun monotoneIncreasing() {
        val g = Grid.geometric(20, 2.0, 5.0, R = 3.0)
        for (i in 0 until g.n) assertTrue(g.x(i + 1) > g.x(i), "node $i not increasing")
    }

    @Test fun hIsMaxStep() {
        val n = 16
        val g = Grid.geometric(n, 0.0, 1.0, R = 2.0)
        var maxStep = 0.0
        for (i in 0 until n) maxStep = maxOf(maxStep, g.x(i + 1) - g.x(i))
        assertEquals(maxStep, g.h, tol)
        // for q>1 the maximal step is the last one.
        assertEquals(g.x(n) - g.x(n - 1), g.h, tol)
    }

    @Test fun extremeStepRatioApproxR() {
        val n = 24; val R = 2.0
        val g = Grid.geometric(n, 0.0, 1.0, R)
        val h0 = g.x(1) - g.x(0)
        val hLast = g.x(n) - g.x(n - 1)
        assertEquals(R, hLast / h0, 1e-9)
    }

    @Test fun localQuasiUniformityRatioEqualsQ() {
        val n = 24; val R = 2.0
        val g = Grid.geometric(n, 0.0, 1.0, R)
        val mu = Math.pow(R, 1.0 / (n - 1)) // = q
        for (i in 1 until n) {
            val ratio = (g.x(i + 1) - g.x(i)) / (g.x(i) - g.x(i - 1))
            assertTrue(ratio >= 1.0 / mu - 1e-9 && ratio <= mu + 1e-9,
                "ratio $ratio out of [mu^-1, mu] at $i")
            assertEquals(mu, ratio, 1e-9)
        }
    }

    @Test fun requiresNAtLeastTwo() {
        assertFailsWith<IllegalArgumentException> { Grid.geometric(1) }
    }
}
