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
 * Граничные случаи библиотеки: минимальные сетки (n = 1, 2), отвержение некорректного входа,
 * поведение вне отрезка сетки и при NaN, пользовательские порождающие системы без локального
 * представления, экстремальные масштабы отрезка и потокобезопасность разделяемых объектов.
 * Зафиксированное поведение экстремальных случаев записывается в `build/reports/boundary-cases.tsv`.
 */
@Tag("fast")
class BoundaryCasesTest {
    private val systems = listOf(GeneratingSystem.B, GeneratingSystem.H, GeneratingSystem.T)

    private companion object {
        /** Файл создаётся заново один раз на запуск JVM, далее только дополняется. */
        private val reportFile: File by lazy {
            File("build/reports/boundary-cases.tsv").also { file ->
                file.parentFile.mkdirs()
                file.writeText("case\toutcome\n")
            }
        }
    }

    private fun report(): File = reportFile

    /** Максимум |sum_j omega_j(t) - 1| по равномерной сетке из m точек отрезка. */
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
            dynamicTest("n=$n, ${sys.name}: ${n + 2} функции, разбиение единицы, точность theta на span phi") {
                val grid = Grid(n, DoubleArray(n + 1) { it.toDouble() / n })
                val basis = MinimalSplineBasis(sys, grid)
                // Ровно n + 2 базисных функции: j = -2..n-1, каждая ненулевая хотя бы в одной точке носителя.
                for (j in -2..n - 1) {
                    val mid = 0.5 * (grid.x(j) + grid.x(j + 3))
                    assertTrue(basis.omega(j, mid) > 0.0, "omega_$j должна быть положительна в середине носителя")
                }
                val pu = partitionOfUnityDefect(basis)
                assertTrue(pu <= 1e-14, "n=$n, ${sys.name}: дефект разбиения единицы $pu > 1e-14")

                val f: (Double) -> Double = { t -> 1.0 + 2.0 * sys.rho(t) + 3.0 * sys.sigma(t) }
                val theta = ProjFunctionals(basis)
                val c = theta.projectorCoeffs(f)
                assertEquals(n + 2, c.size)
                val err = errorEh(f, { t -> basis.evalSpline(c, t) }, grid)
                assertTrue(err <= 1e-12, "n=$n, ${sys.name}: theta не воспроизводит span phi, E_h = $err")
            }
        }
    }

    // ---------------------------------------------------------------- отвержение входа

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
        // NaN нарушает строгое возрастание (любое сравнение с NaN ложно), поэтому отвергается той же проверкой.
        assertFailsWith<IllegalArgumentException> { Grid(2, doubleArrayOf(0.0, Double.NaN, 1.0)) }
        assertFailsWith<IllegalArgumentException> { Grid(1, doubleArrayOf(Double.NaN, 1.0)) }
        assertFailsWith<IllegalArgumentException> { Grid(1, doubleArrayOf(0.0, Double.NaN)) }
    }

    @Test
    fun evalSpline_wrongCoefficientLengthRejected() {
        // Длина вектора коэффициентов проверяется явно: и избыточная, и недостаточная длина
        // отклоняются IllegalArgumentException с указанием ожидаемой длины n + 2.
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
                assertTrue(e.message!!.contains("n + 2 = ${n + 2}") && e.message!!.contains("получено ${c.size}"), e.message)
            }
            // Корректная длина принимается: разбиение единицы даёт 1 для значения и 0 для производных.
            val expected = if (name == "evalSpline") 1.0 else 0.0
            assertEquals(expected, eval(DoubleArray(n + 2) { 1.0 }, 0.5), 1e-12, name)
        }
    }

    // ---------------------------------------------------------------- вне отрезка, NaN

    @Test
    fun omega_indexOutsideRangeRejected() {
        // Индекс сплайна вне [-2, n-1] отклоняется IllegalArgumentException при любом t.
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
                assertTrue(e.message!!.contains("[-2, ${n - 1}]") && e.message!!.contains("получено $j"), e.message)
            }
            // Границы диапазона допустимы, вне носителя значение 0.
            assertEquals(0.0, eval(-2, 1.0))
            assertEquals(0.0, eval(n - 1, 0.0))
        }
        for (k in listOf(-1, n, n + 3)) {
            val e = assertFailsWith<IllegalArgumentException>("activeOmega(k=$k)") { basis.activeOmega(k, 0.5) }
            assertTrue(e.message!!.contains("[0, ${n - 1}]") && e.message!!.contains("получено $k"), e.message)
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
        // Концы отрезка принадлежат ему.
        assertEquals(0, basis.interval(0.0))
        assertEquals(basis.grid.n - 1, basis.interval(1.0))
        // omega вне носителя, в том числе вне отрезка, равна нулю без исключения.
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
            "исключение ${e::class.simpleName}"
        }
        report().appendText("evalSpline(t=NaN)\t$outcome\n")
        // Допустимы оба исхода; сейчас NaN проходит сквозь проверку отрезка и распространяется в результат.
        assertTrue(outcome.startsWith("evalSpline=NaN") || outcome.startsWith("исключение"), outcome)
    }

    // ---------------------------------------------------------------- пользовательская система

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
            dynamicTest("${sys.name} на [${grid.a}, ${grid.b}], n=${grid.n}") {
                val basis = MinimalSplineBasis(sys, grid)
                val pu = partitionOfUnityDefect(basis)
                report().appendText("custom ${sys.name} [${grid.a},${grid.b}] n=${grid.n}\tpu=$pu\n")
                assertTrue(pu <= 1e-10, "${sys.name}: дефект разбиения единицы $pu > 1e-10")
                val f: (Double) -> Double = { t -> 1.0 + 2.0 * sys.rho(t) + 3.0 * sys.sigma(t) }
                val c = ProjFunctionals(basis).projectorCoeffs(f)
                val err = errorEh(f, { t -> basis.evalSpline(c, t) }, grid)
                assertTrue(err <= 1e-9, "${sys.name}: theta не воспроизводит span phi, E_h = $err")
            }
        }
    }

    // ---------------------------------------------------------------- экстремальные отрезки

    @TestFactory
    fun extremeSegments(): List<DynamicTest> = systems.flatMap { sys ->
        listOf(1e6, 1e-6).map { b ->
            dynamicTest("${sys.name} на [0, $b], n=100") {
                val grid = Grid.uniform(100, 0.0, b)
                // H на [0, 1e6] при h = 1e4 не представима в double: локальный масштаб l = min(h, 1) = 1
                // и sinh(2h) переполняется. Ожидается внятное исключение о нефинитности, а не NaN.
                val overflow = sys === GeneratingSystem.H && b > 1.0
                val outcome = try {
                    val basis = MinimalSplineBasis(sys, grid)
                    val pu = partitionOfUnityDefect(basis)
                    val minOmega = (0..400).minOf { i ->
                        val t = b * i / 400
                        (-2..grid.n - 1).minOf { j -> basis.omega(j, t) }
                    }
                    "построен, pu=$pu, min omega=$minOmega"
                } catch (e: IllegalArgumentException) {
                    "IllegalArgumentException: ${e.message?.lines()?.first()?.take(200)}"
                }
                report().appendText("${sys.name} [0,$b] n=100\t$outcome\n")
                if (overflow) {
                    assertTrue(outcome.startsWith("IllegalArgumentException"), "${sys.name} [0,$b]: $outcome")
                    assertTrue(
                        outcome.contains("переполняется") && outcome.contains("нефинитны") && !outcome.contains("NaN"),
                        "${sys.name} [0,$b]: $outcome",
                    )
                } else {
                    assertTrue(outcome.startsWith("построен"), "${sys.name} [0,$b]: $outcome")
                    val pu = outcome.substringAfter("pu=").substringBefore(",").toDouble()
                    assertTrue(pu <= 1e-12, "${sys.name} [0,$b]: дефект разбиения единицы $pu > 1e-12")
                }
            }
        }
    }

    // ---------------------------------------------------------------- потокобезопасность

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
                    assertEquals(refValues[i].toRawBits(), values[i].toRawBits(), "evalSpline в точке ${points[i]} отличается побитово")
                }
            }
            for (fut in coeffFutures) {
                for (coeffs in fut.get()) {
                    assertEquals(refCoeffs.size, coeffs.size)
                    for (i in coeffs.indices) {
                        if (nativeBackend) {
                            val scale = max(1.0, abs(refCoeffs[i]))
                            assertTrue(abs(coeffs[i] - refCoeffs[i]) <= 1e-14 * scale, "projectorCoeffs[$i] расходится: ${coeffs[i]} vs ${refCoeffs[i]}")
                        } else {
                            assertEquals(refCoeffs[i].toRawBits(), coeffs[i].toRawBits(), "projectorCoeffs[$i] отличается побитово")
                        }
                    }
                }
            }
        } finally {
            pool.shutdown()
        }
        report().appendText("threads 8x1000 evalSpline, 8x50 projectorCoeffs\tсовпадают (backend=${Backends.default().name})\n")
    }
}
