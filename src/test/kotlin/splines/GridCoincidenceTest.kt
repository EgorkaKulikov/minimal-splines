package splines

import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Equivalence of the index predicate [Grid.isCoincident] and the value comparison `x(j) == x(j+1)`
 * over the whole admissible index range `j = -2..n+1`, for all four grid factories, for several n
 * and on intervals of different scale and sign; `MinimalSplineBasis.computeA` uses this predicate
 * to detect a multiple boundary node.
 */
@Tag("fast")
class GridCoincidenceTest {

    /** Intervals on which the grids are built (we check independence of scale and sign). */
    private val segments = listOf(0.0 to 1.0, -2.0 to 3.5)

    /**
     * All grids for a given n: uniform and quasiUniform are defined for n >= 1,
     * geometric and graded require n >= 2.
     */
    private fun gridsFor(n: Int, a: Double, b: Double): List<Pair<String, Grid>> = buildList {
        add("uniform" to Grid.uniform(n, a, b))
        add("quasiUniform" to Grid.quasiUniform(n, a, b))
        if (n >= 2) {
            add("geometric" to Grid.geometric(n, a, b))
            add("graded" to Grid.graded(n, a, b))
        }
    }

    /**
     * Main check: `isCoincident(j)` == (`x(j) == x(j+1)`) for all admissible j.
     *
     * The range j = -2..n+1 is the widest one for which both nodes of the pair lie within the
     * stored range -2..n+2. Both ends of the range are included.
     */
    @Test fun predicateMatchesValueComparisonEverywhere() {
        var checks = 0
        var gridsChecked = 0
        var coincidentSeen = 0
        for (n in listOf(1, 2, 8, 16)) {
            for ((a, b) in segments) {
                for ((name, g) in gridsFor(n, a, b)) {
                    gridsChecked++
                    for (j in -2..n + 1) {
                        val byValue = g.x(j) == g.x(j + 1)
                        val byIndex = g.isCoincident(j)
                        assertEquals(
                            byValue, byIndex,
                            "$name(n=$n, [$a,$b]): j=$j -> x($j)=${g.x(j)}, x(${j + 1})=${g.x(j + 1)}; " +
                                "value comparison gives $byValue, the index predicate gives $byIndex",
                        )
                        if (byIndex) coincidentSeen++
                        checks++
                    }
                }
            }
        }
        // Self-check of the test: both branches of the predicate were actually exercised.
        assertEquals(28, gridsChecked, "expected 28 grids (2+4+4+4 factories on 2 intervals)")
        assertEquals(324, checks, "expected 324 checks (n+4 indices per grid)")
        assertEquals(4 * gridsChecked, coincidentSeen, "every grid has exactly 4 coincident pairs: j=-2,-1,n,n+1")
    }

    /**
     * The index shift in `computeA` is checked on `computeA` itself: `computeA(j)` refers to the
     * pair (x_{j+1}, x_{j+2}), i.e. it uses `isCoincident(j + 1)`. The observable consequence of the
     * branching: at a multiple node exactly `phi(x_{j+1})` is returned (with no correction
     * subtracted), while on the general branch the result differs from it (`coef * phiD` is
     * subtracted, `coef != 0`).
     */
    @Test fun computeAUsesCoincidenceAtShiftedIndex() {
        var coincidentSeen = 0
        var regularSeen = 0
        for (n in listOf(2, 8, 16)) {
            for ((a, b) in segments) {
                for ((name, g) in gridsFor(n, a, b)) {
                    val basis = MinimalSplineBasis(GeneratingSystem.B, g)
                    for (j in -2..n - 1) {
                        val actual = basis.computeA(j)
                        val phiAtLeft = basis.sys.phi(g.x(j + 1))
                        val isTripleKnot = g.isCoincident(j + 1)
                        // Self-check of the criterion: the predicate must agree with the old comparison.
                        assertEquals(
                            g.x(j + 1) == g.x(j + 2), isTripleKnot,
                            "$name(n=$n, [$a,$b]): the predicate disagreed with the old criterion at j=$j",
                        )
                        if (isTripleKnot) {
                            coincidentSeen++
                            assertTrue(
                                actual.contentEquals(phiAtLeft),
                                "$name(n=$n, [$a,$b]): at a multiple node (j=$j) computeA must return " +
                                    "exactly phi(x_${j + 1})=${phiAtLeft.toList()}, got ${actual.toList()}",
                            )
                        } else {
                            regularSeen++
                            assertTrue(
                                !actual.contentEquals(phiAtLeft),
                                "$name(n=$n, [$a,$b]): at a NON-multiple node (j=$j) computeA returned exactly " +
                                    "phi(x_${j + 1}) — so the triple-node branch fired by mistake " +
                                    "(likely cause: a wrong index shift in isCoincident)",
                            )
                        }
                    }
                }
            }
        }
        // Guard against a vacuous test: both branches of computeA were actually taken.
        assertTrue(coincidentSeen > 0, "The triple-node branch never fired")
        assertTrue(regularSeen > 0, "The general branch of computeA never fired")
    }

    /**
     * The shift in terms of the grid itself: over the working range of `computeA` (j = -2..n-1)
     * the call `isCoincident(j + 1)` always stays inside the admissible range -2..n+1
     * and agrees with the old value comparison (i.e. the exception is unreachable).
     */
    @Test fun shiftedIndexStaysInAllowedRange() {
        for (n in listOf(1, 2, 8, 16)) {
            for ((a, b) in segments) {
                for ((name, g) in gridsFor(n, a, b)) {
                    for (j in -2..n - 1) {
                        assertTrue(
                            (j + 1) in -2..g.n + 1,
                            "$name(n=$n): the shifted index ${j + 1} left the admissible range",
                        )
                        assertEquals(
                            g.x(j + 1) == g.x(j + 2), g.isCoincident(j + 1),
                            "$name(n=$n, [$a,$b]): the shift is wrong at j=$j",
                        )
                    }
                }
            }
        }
    }

    /** Exactly the boundary pairs are coincident: j = -2, -1 (left triple node) and j = n, n+1 (right). */
    @Test fun onlyBoundaryPairsAreCoincident() {
        val g = Grid.graded(8, 0.0, 1.0)
        for (j in listOf(-2, -1, g.n, g.n + 1)) assertTrue(g.isCoincident(j), "j=$j must be coincident")
        for (j in 0 until g.n) assertTrue(!g.isCoincident(j), "interior j=$j cannot be coincident")
    }

    /** Outside the range -2..n+1 the second node of the pair does not exist — the predicate must reject the input. */
    @Test fun rejectsOutOfRangeIndex() {
        val g = Grid.uniform(4, 0.0, 1.0)
        assertFailsWith<IllegalArgumentException> { g.isCoincident(-3) }
        assertFailsWith<IllegalArgumentException> { g.isCoincident(g.n + 2) }
    }
}
