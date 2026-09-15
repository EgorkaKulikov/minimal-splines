package splines

import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Additional tests of the grid and of the reference formulas: quasiUniform, the require checks,
 * node degeneracy and all branches of ReferenceSplines (B, B', H) inside and outside the support.
 */
@Tag("fast")
class GridExtraTest {
    private val tol = 1e-12

    /** quasiUniform: a monotone grid with the ends a, b and a fixed amplitude. */
    @Test fun quasiUniformMonotoneEnds() {
        val g = Grid.quasiUniform(10, 0.0, 1.0, amp = 0.04)
        assertEquals(0.0, g.a, tol)
        assertEquals(1.0, g.b, tol)
        // strict monotonicity of the interior nodes
        for (i in 0 until g.n) assertTrue(g.x(i + 1) > g.x(i), "node $i not increasing")
        // the boundary is a triple knot
        assertEquals(g.a, g.x(0), tol); assertEquals(g.b, g.x(g.n), tol)
    }

    /** The Grid constructor requires interior of size n+1. */
    @Test fun gridRequiresCorrectSize() {
        assertFailsWith<IllegalArgumentException> { Grid(4, doubleArrayOf(0.0, 1.0)) }
    }

    /**
     * Non-monotone nodes are rejected: without this, h = max_j(x_{j+1}-x_j) and the binary search
     * for the interval would produce unreliable values instead of an error.
     */
    @Test fun gridRejectsNonMonotoneNodes() {
        val e = assertFailsWith<IllegalArgumentException> {
            Grid(3, doubleArrayOf(0.0, 0.7, 0.3, 1.0)) // nodes 1 and 2 are swapped
        }
        assertTrue(e.message!!.contains("i=1"), "the message must point at the index: ${e.message}")
    }

    /** A duplicated interior node (zero step) also violates STRICT increase. */
    @Test fun gridRejectsDuplicatedNode() {
        val e = assertFailsWith<IllegalArgumentException> {
            Grid(3, doubleArrayOf(0.0, 0.5, 0.5, 1.0))
        }
        assertTrue(e.message!!.contains("i=1"), "the message must point at the index: ${e.message}")
    }

    /** n = 0: there are no intervals, the step h = max_j (x_{j+1} - x_j) is undefined; rejected with a diagnostic. */
    @Test fun gridRejectsZeroIntervals() {
        val e = assertFailsWith<IllegalArgumentException> { Grid(0, doubleArrayOf(0.0)) }
        assertTrue(e.message!!.contains("n=0"), "the message must contain n: ${e.message}")
    }

    /**
     * All four factories build STRICTLY increasing nodes for every n and every parameter in use —
     * that is, the new monotonicity check does not fire on the regular paths.
     */
    @Test fun allFactoriesProduceStrictlyIncreasingNodes() {
        val ns = listOf(2, 3, 4, 5, 8, 16, 32, 64, 128)
        for (n in ns) {
            val grids = buildList {
                add("uniform(n=$n)" to Grid.uniform(n, 0.0, 1.0))
                add("uniform(n=$n,[-3,7])" to Grid.uniform(n, -3.0, 7.0))
                for (amp in listOf(-0.04, 0.0, 0.04, 0.15)) {
                    add("quasiUniform(n=$n,amp=$amp)" to Grid.quasiUniform(n, 0.0, 1.0, amp))
                }
                for (r in listOf(0.1, 0.5, 1.0, 1.5, 2.0, 10.0, 1e6)) {
                    add("geometric(n=$n,R=$r)" to Grid.geometric(n, 0.0, 1.0, r))
                    add("graded(n=$n,ratio=$r)" to Grid.graded(n, 0.0, 1.0, r))
                }
            }
            for ((name, g) in grids) {
                for (i in 0 until g.n) {
                    assertTrue(g.x(i + 1) > g.x(i), "$name: nodes do not increase at i=$i")
                }
                assertTrue(g.h > 0.0, "$name: h <= 0")
            }
        }
    }

    /**
     * quasiUniform is the only factory able to produce non-monotone nodes:
     * Psi'(u) = 1 + 2*pi*amp*cos(2*pi*u) changes sign for |amp| > 1/(2*pi) ≈ 0.15915.
     *
     * The rejection is driven by the ACTUAL nodes (the [Grid] invariant) and NOT by amp itself:
     * non-monotonicity of Psi as a function is not equivalent to non-monotonicity of a finite set
     * of nodes. Both directions are checked: non-monotone inputs are rejected, while admissible ones
     * (amp = 1/(2*pi) for any n; amp = 0.16 for n = 8) are built.
     */
    @Test fun quasiUniformRejectsOnlyActuallyNonMonotoneNodes() {
        val e = assertFailsWith<IllegalArgumentException> { Grid.quasiUniform(19, 0.0, 1.0, amp = 0.16) }
        assertTrue(
            e.message!!.contains("strictly increasing"),
            "the message must point at the non-monotone nodes: ${e.message}",
        )
        assertFailsWith<IllegalArgumentException> { Grid.quasiUniform(8, 0.0, 1.0, amp = -0.5) }
        // The same amp = 0.16 with n = 8 gives strictly increasing nodes and is admissible.
        val coarse = Grid.quasiUniform(8, 0.0, 1.0, amp = 0.16)
        for (i in 0 until coarse.n) assertTrue(coarse.x(i + 1) > coarse.x(i), "amp=0.16,n=8: node $i")
        // Exactly at the boundary 1/(2*pi): Psi' vanishes at a SINGLE point, but the nodes are
        // strictly increasing for any finite n — there is nothing to reject.
        for (n in intArrayOf(8, 64, 1024)) {
            val g = Grid.quasiUniform(n, 0.0, 1.0, amp = 1.0 / (2.0 * Math.PI))
            for (i in 0 until g.n) assertTrue(g.x(i + 1) > g.x(i), "amp=1/(2pi),n=$n: node $i")
        }
        // slightly below the boundary — still built
        val g = Grid.quasiUniform(64, 0.0, 1.0, amp = 0.158)
        for (i in 0 until g.n) assertTrue(g.x(i + 1) > g.x(i), "amp=0.158: node $i")
    }

    /** nonDegenerate=false when nodes coincide at the triple knot on the left edge (j=-2). */
    @Test fun nonDegenerateFalseAtTripleKnot() {
        val g = Grid.uniform(8)
        // x(-2)=x(-1)=x(0)=a -> degenerate
        assertTrue(!nonDegenerate(g, -2))
    }

    /** omegaB is zero outside the support [x_j, x_{j+3}] (both branches t<xj and t>xj3). */
    @Test fun omegaBZeroOutsideSupport() {
        val g = Grid.uniform(8)
        val j = 2
        assertEquals(0.0, ReferenceSplines.omegaB(g, j, g.x(j) - 0.01), tol)
        assertEquals(0.0, ReferenceSplines.omegaB(g, j, g.x(j + 3) + 0.01), tol)
        assertEquals(0.0, ReferenceSplines.omegaBDeriv(g, j, g.x(j) - 0.01), tol)
        assertEquals(0.0, ReferenceSplines.omegaBDeriv(g, j, g.x(j + 3) + 0.01), tol)
    }

    /** omegaB is non-negative and covers all three pieces (left, middle, right). */
    @Test fun omegaBThreePieces() {
        val g = Grid.uniform(8)
        val j = 2
        val tLeft = 0.5 * (g.x(j) + g.x(j + 1))
        val tMid = 0.5 * (g.x(j + 1) + g.x(j + 2))
        val tRight = 0.5 * (g.x(j + 2) + g.x(j + 3))
        for (t in listOf(tLeft, tMid, tRight)) {
            assertTrue(ReferenceSplines.omegaB(g, j, t) > 0.0, "omegaB<=0 at $t")
        }
    }

    /** omegaBDeriv — the numerical derivative agrees with the analytic one on every piece. */
    @Test fun omegaBDerivMatchesNumeric() {
        val g = Grid.uniform(8)
        val j = 2
        val eps = 1e-6
        for (t in listOf(g.x(j) + 0.02, 0.5 * (g.x(j + 1) + g.x(j + 2)), g.x(j + 3) - 0.02)) {
            val num = (ReferenceSplines.omegaB(g, j, t + eps) - ReferenceSplines.omegaB(g, j, t - eps)) / (2 * eps)
            assertEquals(num, ReferenceSplines.omegaBDeriv(g, j, t), 1e-4, "deriv mismatch at $t")
        }
    }

    /** omegaH: zero outside the support and positive on all three pieces. */
    @Test fun omegaHPiecesAndSupport() {
        val g = Grid.uniform(8)
        val j = 2
        assertEquals(0.0, ReferenceSplines.omegaH(g, j, g.x(j) - 0.01), tol)
        assertEquals(0.0, ReferenceSplines.omegaH(g, j, g.x(j + 3) + 0.01), tol)
        val tLeft = 0.5 * (g.x(j) + g.x(j + 1))
        val tMid = 0.5 * (g.x(j + 1) + g.x(j + 2))
        val tRight = 0.5 * (g.x(j + 2) + g.x(j + 3))
        for (t in listOf(tLeft, tMid, tRight)) {
            assertTrue(ReferenceSplines.omegaH(g, j, t) > 0.0, "omegaH<=0 at $t")
        }
    }
}
