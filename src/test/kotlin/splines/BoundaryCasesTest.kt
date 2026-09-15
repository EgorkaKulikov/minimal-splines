package splines

import numerics.backend.Backends
import org.junit.jupiter.api.DynamicTest
import org.junit.jupiter.api.DynamicTest.dynamicTest
import org.junit.jupiter.api.Tag
import org.junit.jupiter.api.TestFactory
import splines.functionals.ProjFunctionals
import splines.metrics.errorEh
import java.io.File
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import kotlin.math.abs
import kotlin.math.exp
import kotlin.math.max
import kotlin.math.sin
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Boundary cases of the library: minimal grids (n = 1, 2), rejection of invalid input, behaviour
 * outside the grid interval and at NaN, user-defined generating systems without a local
 * representation, extreme interval scales and thread safety of shared objects.
 * The recorded behaviour of the extreme cases is written to `build/reports/boundary-cases.tsv`.
 */
@Tag("fast")
class BoundaryCasesTest {
    private val systems = listOf(GeneratingSystem.B, GeneratingSystem.H, GeneratingSystem.T)

    private companion object {
        /** The file is created from scratch once per JVM run and is only appended to afterwards. */
        private val reportFile: File by lazy {
            File("build/reports/boundary-cases.tsv").also { file ->
                file.parentFile.mkdirs()
                file.writeText("case\toutcome\n")
            }
        }
    }

    private fun report(): File = reportFile

    /** Maximum of |sum_j omega_j(t) - 1| over a uniform grid of m points of the interval. */
    private fun partitionOfUnityDefect(basis: MinimalSplineBasis, m: Int = 400): Double {
        val g = basis.grid
        var worst = 0.0
        for (i in 0..m) {
            val t = g.a + (g.b - g.a) * i / m
            var s = 0.0
            for (j in -2..g.n - 1) s += basis.omega(j, t)
            worst = max(worst, abs(s - 1.0))
        }
        return worst
    }

    // ---------------------------------------------------------------- n = 1, n = 2

    @TestFactory
    fun smallestGrids(): List<DynamicTest> = listOf(1, 2).flatMap { n ->
        systems.map { sys ->
            dynamicTest("n=$n, ${sys.name}: ${n + 2} functions, partition of unity, accuracy of theta on span phi") {
                val grid = Grid(n, DoubleArray(n + 1) { it.toDouble() / n })
                val basis = MinimalSplineBasis(sys, grid)
                // Exactly n + 2 basis functions: j = -2..n-1, each nonzero at least at one point of its support.
                for (j in -2..n - 1) {
                    val mid = 0.5 * (grid.x(j) + grid.x(j + 3))
                    assertTrue(basis.omega(j, mid) > 0.0, "omega_$j must be positive at the middle of its support")
                }
                val pu = partitionOfUnityDefect(basis)
                assertTrue(pu <= 1e-14, "n=$n, ${sys.name}: partition of unity defect $pu > 1e-14")

                val f: (Double) -> Double = { t -> 1.0 + 2.0 * sys.rho(t) + 3.0 * sys.sigma(t) }
                val theta = ProjFunctionals(basis)
                val c = theta.projectorCoeffs(f)
                assertEquals(n + 2, c.size)
                val err = errorEh(f, { t -> basis.evalSpline(c, t) }, grid)
                assertTrue(err <= 1e-12, "n=$n, ${sys.name}: theta does not reproduce span phi, E_h = $err")
            }
        }
    }

    // ---------------------------------------------------------------- input rejection

    @Test
    fun grid_nonIncreasingInteriorRejected() {
        assertFailsWith<IllegalArgumentException> { Grid(3, doubleArrayOf(0.0, 0.5, 0.5, 1.0)) }
        assertFailsWith<IllegalArgumentException> { Grid(3, doubleArrayOf(0.0, 0.7, 0.3, 1.0)) }
    }

    @Test
    fun grid_sizeMismatchRejected() {
        assertFailsWith<IllegalArgumentException> { Grid(3, doubleArrayOf(0.0, 0.5, 1.0)) }
        assertFailsWith<IllegalArgumentException> { Grid(2, doubleArrayOf(0.0, 0.25, 0.5, 1.0)) }
    }

    @Test
    fun grid_reversedOrEmptyIntervalRejected() {
        assertFailsWith<IllegalArgumentException> { Grid.uniform(8, 1.0, 0.0) }
        assertFailsWith<IllegalArgumentException> { Grid.uniform(8, 0.0, 0.0) }
        assertFailsWith<IllegalArgumentException> { Grid.graded(8, 1.0, 0.0) }
    }

    @Test
    fun grid_zeroIntervalsRejected() {
        assertFailsWith<IllegalArgumentException> { Grid(0, doubleArrayOf(0.0)) }
        assertFailsWith<IllegalArgumentException> { Grid.uniform(0, 0.0, 1.0) }
    }

    @Test
    fun grid_nanInteriorRejected() {
        // NaN violates strict increase (any comparison with NaN is false), so it is rejected by the same check.
        assertFailsWith<IllegalArgumentException> { Grid(2, doubleArrayOf(0.0, Double.NaN, 1.0)) }
        assertFailsWith<IllegalArgumentException> { Grid(1, doubleArrayOf(Double.NaN, 1.0)) }
        assertFailsWith<IllegalArgumentException> { Grid(1, doubleArrayOf(0.0, Double.NaN)) }
    }

    @Test
    fun evalSpline_wrongCoefficientLengthRejected() {
        // The length of the coefficient vector is checked explicitly: both an excessive and an
        // insufficient length are rejected with IllegalArgumentException stating the expected length n + 2.
        val basis = MinimalSplineBasis(GeneratingSystem.B, Grid.uniform(8, 0.0, 1.0))
        val n = basis.grid.n
        val long = DoubleArray(n + 3) { 1.0 }
        val short = DoubleArray(n + 1) { 1.0 }
        val evaluators = listOf<Pair<String, (DoubleArray, Double) -> Double>>(
            "evalSpline" to basis::evalSpline,
            "evalSplineDeriv" to basis::evalSplineDeriv,
            "evalSplineDeriv2" to basis::evalSplineDeriv2,
        )
        for ((name, eval) in evaluators) {
            for ((c, t) in listOf(long to 0.5, short to 0.5, short to 0.05, short to 1.0, DoubleArray(0) to 0.5)) {
                val e = assertFailsWith<IllegalArgumentException>("$name(c.size=${c.size}, t=$t)") { eval(c, t) }
                assertTrue(e.message!!.contains("n + 2 = ${n + 2}") && e.message!!.contains("got ${c.size}"), e.message)
            }
            // A correct length is accepted: the partition of unity gives 1 for the value and 0 for the derivatives.
            val expected = if (name == "evalSpline") 1.0 else 0.0
            assertEquals(expected, eval(DoubleArray(n + 2) { 1.0 }, 0.5), 1e-12, name)
        }
    }

    // ---------------------------------------------------------------- outside the interval, NaN

    @Test
    fun omega_indexOutsideRangeRejected() {
        // A spline index outside [-2, n-1] is rejected with IllegalArgumentException for any t.
        val basis = MinimalSplineBasis(GeneratingSystem.B, Grid.uniform(8, 0.0, 1.0))
        val n = basis.grid.n
        val evaluators = listOf<Pair<String, (Int, Double) -> Double>>(
            "omega" to basis::omega,
            "omegaDeriv" to basis::omegaDeriv,
            "omegaDeriv2" to basis::omegaDeriv2,
        )
        for ((name, eval) in evaluators) {
            for ((j, t) in listOf(-3 to 0.5, n to 0.5, n to 1.0, n + 2 to 0.5, n + 5 to 0.5)) {
                val e = assertFailsWith<IllegalArgumentException>("$name(j=$j, t=$t)") { eval(j, t) }
                assertTrue(e.message!!.contains("[-2, ${n - 1}]") && e.message!!.contains("got $j"), e.message)
            }
            // The ends of the range are admissible; outside the support the value is 0.
            assertEquals(0.0, eval(-2, 1.0))
            assertEquals(0.0, eval(n - 1, 0.0))
        }
        for (k in listOf(-1, n, n + 3)) {
            val e = assertFailsWith<IllegalArgumentException>("activeOmega(k=$k)") { basis.activeOmega(k, 0.5) }
            assertTrue(e.message!!.contains("[0, ${n - 1}]") && e.message!!.contains("got $k"), e.message)
        }
        assertEquals(3, basis.activeOmega(0, 0.0).size)
        assertEquals(3, basis.activeOmega(n - 1, 1.0).size)
    }

    @Test
    fun evaluationOutsideSegmentRejected_noExtrapolation() {
        val basis = MinimalSplineBasis(GeneratingSystem.T, Grid.uniform(8, 0.0, 1.0))
        val c = DoubleArray(basis.grid.n + 2) { 1.0 + it }
        for (t in doubleArrayOf(-1e-9, -0.5, 1.0 + 1e-9, 2.0)) {
            assertFailsWith<IllegalArgumentException>("t=$t") { basis.interval(t) }
            assertFailsWith<IllegalArgumentException>("t=$t") { basis.evalSpline(c, t) }
            assertFailsWith<IllegalArgumentException>("t=$t") { basis.evalSplineDeriv(c, t) }
            assertFailsWith<IllegalArgumentException>("t=$t") { basis.evalSplineDeriv2(c, t) }
        }
        // The ends of the interval belong to it.
        assertEquals(0, basis.interval(0.0))
        assertEquals(basis.grid.n - 1, basis.interval(1.0))
        // omega outside its support, including outside the interval, is zero without an exception.
        assertEquals(0.0, basis.omega(0, -1.0))
        assertEquals(0.0, basis.omega(basis.grid.n - 1, 2.0))
    }

    @Test
    fun evaluationAtNan_recordCurrentBehaviour() {
        val basis = MinimalSplineBasis(GeneratingSystem.B, Grid.uniform(8, 0.0, 1.0))
        val c = DoubleArray(basis.grid.n + 2) { 1.0 }
        val outcome = try {
            val v = basis.evalSpline(c, Double.NaN)
            val w = basis.omega(0, Double.NaN)
            "evalSpline=$v, omega=$w"
        } catch (e: RuntimeException) {
            "exception ${e::class.simpleName}"
        }
        report().appendText("evalSpline(t=NaN)\t$outcome\n")
        // Both outcomes are admissible; currently NaN passes through the interval check and propagates into the result.
        assertTrue(outcome.startsWith("evalSpline=NaN") || outcome.startsWith("exception"), outcome)
    }

    // ---------------------------------------------------------------- user-defined system

    @TestFactory
    fun customSystemWithoutLocalFrame(): List<DynamicTest> {
        val expSys = GeneratingSystem(
            "exp", rho = { t -> exp(t) }, sigma = { t -> exp(2.0 * t) },
            rhoD = { t -> exp(t) }, sigmaD = { t -> 2.0 * exp(2.0 * t) },
            rhoDD = { t -> exp(t) }, sigmaDD = { t -> 4.0 * exp(2.0 * t) },
        )
        val cubicSys = GeneratingSystem(
            "cubic", rho = { t -> t }, sigma = { t -> t * t * t },
            rhoD = { 1.0 }, sigmaD = { t -> 3.0 * t * t },
            rhoDD = { 0.0 }, sigmaDD = { t -> 6.0 * t },
        )
        return listOf(expSys to Grid.uniform(50, 0.0, 1.0), cubicSys to Grid.uniform(50, 1.0, 2.0)).map { (sys, grid) ->
            dynamicTest("${sys.name} on [${grid.a}, ${grid.b}], n=${grid.n}") {
                val basis = MinimalSplineBasis(sys, grid)
                val pu = partitionOfUnityDefect(basis)
                report().appendText("custom ${sys.name} [${grid.a},${grid.b}] n=${grid.n}\tpu=$pu\n")
                assertTrue(pu <= 1e-10, "${sys.name}: partition of unity defect $pu > 1e-10")
                val f: (Double) -> Double = { t -> 1.0 + 2.0 * sys.rho(t) + 3.0 * sys.sigma(t) }
                val c = ProjFunctionals(basis).projectorCoeffs(f)
                val err = errorEh(f, { t -> basis.evalSpline(c, t) }, grid)
                assertTrue(err <= 1e-9, "${sys.name}: theta does not reproduce span phi, E_h = $err")
            }
        }
    }

    // ---------------------------------------------------------------- extreme intervals

    @TestFactory
    fun extremeSegments(): List<DynamicTest> = systems.flatMap { sys ->
        listOf(1e6, 1e-6).map { b ->
            dynamicTest("${sys.name} on [0, $b], n=100") {
                val grid = Grid.uniform(100, 0.0, b)
                // H on [0, 1e6] with h = 1e4 is not representable in double: the local scale is l = min(h, 1) = 1
                // and sinh(2h) overflows. A clear exception about non-finite values is expected, not a NaN.
                val overflow = sys === GeneratingSystem.H && b > 1.0
                val outcome = try {
                    val basis = MinimalSplineBasis(sys, grid)
                    val pu = partitionOfUnityDefect(basis)
                    val minOmega = (0..400).minOf { i ->
                        val t = b * i / 400
                        (-2..grid.n - 1).minOf { j -> basis.omega(j, t) }
                    }
                    "built, pu=$pu, min omega=$minOmega"
                } catch (e: IllegalArgumentException) {
                    "IllegalArgumentException: ${e.message?.lines()?.first()?.take(200)}"
                }
                report().appendText("${sys.name} [0,$b] n=100\t$outcome\n")
                if (overflow) {
                    assertTrue(outcome.startsWith("IllegalArgumentException"), "${sys.name} [0,$b]: $outcome")
                    assertTrue(
                        outcome.contains("overflows on interval") && outcome.contains("are non-finite") && !outcome.contains("NaN"),
                        "${sys.name} [0,$b]: $outcome",
                    )
                } else {
                    assertTrue(outcome.startsWith("built"), "${sys.name} [0,$b]: $outcome")
                    val pu = outcome.substringAfter("pu=").substringBefore(",").toDouble()
                    assertTrue(pu <= 1e-12, "${sys.name} [0,$b]: partition of unity defect $pu > 1e-12")
                }
            }
        }
    }

    // ---------------------------------------------------------------- thread safety

    @Test
    fun sharedBasisAndFunctionals_threadSafe() {
        val basis = MinimalSplineBasis(GeneratingSystem.B, Grid.uniform(64, 0.0, 1.0))
        val theta = ProjFunctionals(basis)
        val n = basis.grid.n
        val c = DoubleArray(n + 2) { sin(1.7 * it) }
        val points = DoubleArray(1000) { (it + 0.5) / 1000.0 }
        val refValues = DoubleArray(points.size) { basis.evalSpline(c, points[it]) }
        val f: (Double) -> Double = { t -> exp(sin(3.0 * t)) }
        val refCoeffs = theta.projectorCoeffs(f)
        val nativeBackend = theta.ctx.backend.isNative

        val pool = Executors.newFixedThreadPool(8)
        try {
            val evalFutures = (0 until 8).map {
                pool.submit(Callable { DoubleArray(points.size) { i -> basis.evalSpline(c, points[i]) } })
            }
            val coeffFutures = (0 until 8).map {
                pool.submit(Callable { List(50) { theta.projectorCoeffs(f) } })
            }
            for (fut in evalFutures) {
                val values = fut.get()
                for (i in points.indices) {
                    assertEquals(refValues[i].toRawBits(), values[i].toRawBits(), "evalSpline at the point ${points[i]} differs bitwise")
                }
            }
            for (fut in coeffFutures) {
                for (coeffs in fut.get()) {
                    assertEquals(refCoeffs.size, coeffs.size)
                    for (i in coeffs.indices) {
                        if (nativeBackend) {
                            val scale = max(1.0, abs(refCoeffs[i]))
                            assertTrue(abs(coeffs[i] - refCoeffs[i]) <= 1e-14 * scale, "projectorCoeffs[$i] diverges: ${coeffs[i]} vs ${refCoeffs[i]}")
                        } else {
                            assertEquals(refCoeffs[i].toRawBits(), coeffs[i].toRawBits(), "projectorCoeffs[$i] differs bitwise")
                        }
                    }
                }
            }
        } finally {
            pool.shutdown()
        }
        report().appendText("threads 8x1000 evalSpline, 8x50 projectorCoeffs\tmatch (backend=${Backends.default().name})\n")
    }
}
