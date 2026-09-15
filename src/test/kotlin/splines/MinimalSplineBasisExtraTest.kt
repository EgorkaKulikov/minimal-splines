package splines

import org.junit.jupiter.api.Tag
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import kotlin.test.fail

/**
 * Additional tests of the minimal spline basis: interval/activeOmega, evaluation of the spline and
 * of its derivative, agreement of omegaDeriv with the numerical derivative, and node degeneracy at
 * the right triple knot.
 */
@Tag("fast")
class MinimalSplineBasisExtraTest {
    private val grid = Grid.uniform(8)
    private val b = MinimalSplineBasis(GeneratingSystem.B, grid)

    /** interval(t): returns the index k with x_k<=t<x_{k+1}; for t=b it returns n-1. */
    @Test fun intervalIndexing() {
        assertEquals(0, b.interval(grid.x(0) + 1e-6))
        assertEquals(grid.n - 1, b.interval(grid.b))
        val k = b.interval(0.5 * (grid.x(3) + grid.x(4)))
        assertEquals(3, k)
    }

    /** activeOmega: the three active splines sum to 1 (partition of unity for B). */
    @Test fun activeOmegaPartitionOfUnity() {
        val t = 0.43
        val k = b.interval(t)
        val w = b.activeOmega(k, t)
        assertEquals(1.0, w[0] + w[1] + w[2], 1e-10)
    }

    /** evalSpline reproduces the constant 1 (sum c_j omega_j with c_j=1 is the partition of unity). */
    @Test fun evalSplineReproducesConstant() {
        val c = DoubleArray(grid.n + 2) { 1.0 }
        for (t in listOf(0.1, 0.37, 0.62, 0.88)) {
            assertEquals(1.0, b.evalSpline(c, t), 1e-9, "const not reproduced at $t")
        }
    }

    /** evalSpline reproduces rho(t)=t for B: the coefficients are the values of rho at the de Boor nodes. */
    @Test fun evalSplineReproducesLinear() {
        // For phi^B=(1,t,t^2) the spline reproduces a linear function exactly.
        // Take the projector theta, project f(t)=t and check the reproduction.
        val theta = splines.functionals.ProjFunctionals(b)
        val c = theta.projectorCoeffs({ t -> t })
        for (t in listOf(0.15, 0.5, 0.81)) {
            assertEquals(t, b.evalSpline(c, t), 1e-8, "linear not reproduced at $t")
        }
    }

    /** evalSplineDeriv: the derivative of the linear function t equals 1. */
    @Test fun evalSplineDerivOfLinear() {
        val theta = splines.functionals.ProjFunctionals(b)
        val c = theta.projectorCoeffs({ t -> t })
        for (t in listOf(0.2, 0.55, 0.77)) {
            assertEquals(1.0, b.evalSplineDeriv(c, t), 1e-7, "deriv != 1 at $t")
        }
    }

    /** omegaDeriv: agrees with the central numerical derivative of omega. */
    @Test fun omegaDerivMatchesNumeric() {
        val j = 2
        val eps = 1e-6
        for (t in listOf(grid.x(j) + 0.03, 0.5 * (grid.x(j + 1) + grid.x(j + 2)), grid.x(j + 3) - 0.03)) {
            val num = (b.omega(j, t + eps) - b.omega(j, t - eps)) / (2 * eps)
            assertEquals(num, b.omegaDeriv(j, t), 1e-4, "omegaDeriv mismatch at $t")
        }
    }

    /** omega and omegaDeriv are zero outside the support [x_j, x_{j+3}]. */
    @Test fun omegaZeroOutsideSupport() {
        val j = 3
        assertEquals(0.0, b.omega(j, grid.x(j) - 0.01), 1e-12)
        assertEquals(0.0, b.omega(j, grid.x(j + 3) + 0.01), 1e-12)
        assertEquals(0.0, b.omegaDeriv(j, grid.x(j) - 0.01), 1e-12)
        assertEquals(0.0, b.omegaDeriv(j, grid.x(j + 3) + 0.01), 1e-12)
    }

    /** omega returns 0 when slot is outside [0,2] (the point lies in the support, but in another active interval). */
    @Test fun omegaZeroWhenSlotOutOfRange() {
        // omega_{-2} has the support [x_{-2},x_1]=[a, x_1]; on the interval k=0 slot=j-(k-2)= -2+2=0 is fine,
        // but for j whose support lies to the right of the current interval slot<0 -> 0.
        val t = 0.5 * (grid.x(0) + grid.x(1)) // interval k=0, active j=-2,-1,0
        assertEquals(0.0, b.omega(1, t), 1e-12) // is j=1 in the support? no: x_1>t -> already 0 by the support
        // explicit check of slot>2: take a j significantly further to the left
        val t2 = 0.5 * (grid.x(5) + grid.x(6)) // interval k=5, active 3,4,5
        assertTrue(abs(b.omega(2, t2)) < 1e-12) // slot=2-(5-2)=-1 <0
    }

    /** nonDegenerate=false at the right triple knot: the first part is true, the second is false. */
    @Test fun nonDegenerateFalseRightEdge() {
        // j=n-1: x(n-1)<x(n) holds, but x(n)==x(n+1) -> the second comparison is false
        assertTrue(!nonDegenerate(grid, grid.n - 1))
    }

    /**
     * The binary interval search gives the same result as the linear one at every point of the
     * interval (inside the intervals, at the interior nodes, at the ends). Points outside the interval
     * are not included in the set: the reference linear search clamps them to the outermost intervals,
     * whereas the contract of [MinimalSplineBasis.interval] is an exception; it is checked in
     * [outsideSegmentIsRejected].
     */
    @Test fun intervalBinarySearchMatchesLinear() {
        val g = Grid.uniform(5)
        val basis = MinimalSplineBasis(GeneratingSystem.B, g)
        // Reference linear search.
        fun linear(t: Double): Int {
            var k = 0
            while (k < g.n - 1 && t >= g.x(k + 1)) k++
            return k
        }
        val pts = mutableListOf<Double>()
        pts.add(g.a)                     // left end
        pts.add(g.b)                     // right end
        for (i in 0..g.n) {
            pts.add(g.x(i))              // exactly at a node (including the interior ones)
            if (i < g.n) pts.add(0.5 * (g.x(i) + g.x(i + 1))) // midpoint of the interval
        }
        for (t in pts) {
            assertEquals(linear(t), basis.interval(t), "interval mismatch at t=$t")
        }
    }

    /**
     * A point outside the interval is rejected with an exception rather than clamped to the outermost
     * interval: outside the interval the spline is not defined, and extrapolating the polynomial of
     * the outermost layer is inadmissible. All four public entry points are checked, since the check
     * lives in the shared `intervalOf` and the test records that each of them goes through it.
     */
    @Test fun outsideSegmentIsRejected() {
        val g = Grid.uniform(5)
        val basis = MinimalSplineBasis(GeneratingSystem.B, g)
        val c = DoubleArray(g.n + 2) { 1.0 }
        val outside = listOf(
            g.a - 0.1 to "0.1 to the left of a",
            g.b + 0.1 to "0.1 to the right of b",
            Math.nextAfter(g.a, Double.NEGATIVE_INFINITY) to "nextDown(a) — one ulp to the left of a",
            Math.nextAfter(g.b, Double.POSITIVE_INFINITY) to "nextUp(b) — one ulp to the right of b",
        )
        val entryPoints: List<Pair<String, (Double) -> Any>> = listOf(
            "interval" to { t -> basis.interval(t) },
            "evalSpline" to { t -> basis.evalSpline(c, t) },
            "evalSplineDeriv" to { t -> basis.evalSplineDeriv(c, t) },
            "evalSplineDeriv2" to { t -> basis.evalSplineDeriv2(c, t) },
        )
        for ((t, where) in outside) {
            for ((name, call) in entryPoints) {
                val error = try {
                    val value = call(t)
                    fail(
                        "$name(t=$t) ($where, interval [${g.a}, ${g.b}]) must have thrown " +
                            "IllegalArgumentException, but returned $value — that is, it extrapolated " +
                            "beyond the interval",
                    )
                } catch (e: IllegalArgumentException) {
                    e
                }
                val message = error.message ?: ""
                // The message must be diagnostic: it must name the point, both bounds and explain
                // that outside the interval the spline is not defined — otherwise the exception
                // merely replaces one obscure outcome with another.
                for (fragment in listOf("MinimalSplineBasis", "$t", "${g.a}", "${g.b}", "not extrapolated")) {
                    assertTrue(
                        message.contains(fragment),
                        "The error message of $name is uninformative: the fragment \"$fragment\" is missing. " +
                            "Got: \"$message\"",
                    )
                }
            }
        }
    }

    /**
     * `t = NaN` does not lead to an exception: the membership check is written with negations, and NaN
     * propagates into the result. Writing the check as `require(t in a..b)` would turn NaN into an
     * exception.
     */
    @Test fun nanArgumentPassesThroughWithoutThrowing() {
        val g = Grid.uniform(5)
        val c = DoubleArray(g.n + 2) { 1.0 }
        for (system in listOf(GeneratingSystem.B, GeneratingSystem.H, GeneratingSystem.T)) {
            val basis = MinimalSplineBasis(system, g)
            // None of the entry points raises an exception for NaN.
            val k = basis.interval(Double.NaN)
            assertTrue(k in 0 until g.n, "${system.name}: interval(NaN) returned $k outside [0, ${g.n - 1}]")
            val value = basis.evalSpline(c, Double.NaN)
            val deriv = basis.evalSplineDeriv(c, Double.NaN)
            val deriv2 = basis.evalSplineDeriv2(c, Double.NaN)

            // Where phi depends on t, NaN propagates into the result.
            // The value and the first derivative behave this way for all three systems.
            assertTrue(value.isNaN(), "${system.name}: evalSpline(NaN) must give NaN, got $value")
            assertTrue(deriv.isNaN(), "${system.name}: evalSplineDeriv(NaN) must give NaN, got $deriv")

            // The second derivative is a special case. For the polynomial phi^B = (1, t, t^2)
            // the second derivative phi'' = (0, 0, 2) is constant: the argument does not enter the
            // computation, so the result is finite for any t,
            // including NaN. For H (sinh/cosh) and T (sin/cos) phi'' depends on t, and NaN propagates.
            if (system.name == "B") {
                assertTrue(
                    deriv2.isFinite(),
                    "B: phi''=(0,0,2) is constant, therefore evalSplineDeriv2(NaN) must be " +
                        "finite (the argument does not take part in the computation), got $deriv2",
                )
            } else {
                assertTrue(
                    deriv2.isNaN(),
                    "${system.name}: phi'' depends on t, therefore evalSplineDeriv2(NaN) must give " +
                        "NaN, got $deriv2",
                )
            }
        }
    }

    /**
     * The `omega*` methods return 0 without an exception at points outside the interval: they are
     * preceded by a support cut-off, and the support of any basis spline lies inside the interval. The
     * outermost j, where the support touches the boundary (j = -2 at `a`, j = n-1 at `b`), are checked
     * as well.
     */
    @Test fun omegaFamilyOutsideSegmentReturnsZeroWithoutThrowing() {
        val g = Grid.uniform(5)
        val basis = MinimalSplineBasis(GeneratingSystem.B, g)
        val outside = listOf(g.a - 0.1, g.b + 0.1, g.a - 1e-12, g.b + 1e-12)
        for (j in -2 until g.n) {
            for (t in outside) {
                assertEquals(0.0, basis.omega(j, t), 0.0, "omega($j, $t) outside the interval must be 0")
                assertEquals(0.0, basis.omegaDeriv(j, t), 0.0, "omegaDeriv($j, $t) outside the interval must be 0")
                assertEquals(0.0, basis.omegaDeriv2(j, t), 0.0, "omegaDeriv2($j, $t) outside the interval must be 0")
            }
        }
    }
}
