package splines

import numerics.orders
import org.junit.jupiter.api.DynamicTest
import org.junit.jupiter.api.DynamicTest.dynamicTest
import org.junit.jupiter.api.Tag
import org.junit.jupiter.api.TestFactory
import splines.functionals.AveragingFunctionals
import splines.functionals.DeBoorFixFunctionals
import splines.functionals.DiscreteDeBoorFixFunctionals
import splines.functionals.FunctionalFamily
import splines.functionals.ProjFunctionals
import splines.functionals.ThreePointFunctionals
import splines.metrics.DEFAULT_CONTROL_REFINEMENT
import splines.metrics.errorEh
import java.io.File
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.exp
import kotlin.math.max
import kotlin.math.sin
import kotlin.test.assertTrue

/**
 * Порядки сходимости квазиинтерполянтов на квадратичных минимальных сплайнах при измельчении
 * сетки: n = 8, 16, 32, 64, 128 на [0, 1], тестовая функция f(t) = exp(sin 3t), не лежащая
 * в span φ ни одной из систем B, H, T.
 *
 * Проверяется:
 *  - монотонное убывание E_h;
 *  - порядок по значению на парах 32→64 и 64→128 в окне [2.7, 3.3] для θ, ξ, μ, λ
 *    (теоретический порядок 3);
 *  - для ξ̃ — порядок 3 во внутренней области (без краевого слоя в 3h) и порядок 2 по полной
 *    норме: односторонняя разность у кратных краевых узлов даёт краевой дефект O(h²);
 *  - порядок по производной для θ и ξ в окне [1.7, 2.3] (теоретический порядок 2).
 *
 * Все измерения записываются в `build/reports/convergence-orders.tsv`
 * (колонки `sys family grid quantity n Eh order`).
 */
@Tag("fast")
class ConvergenceOrderTest {
    private val ns = listOf(8, 16, 32, 64, 128)

    private val f: (Double) -> Double = { t -> exp(sin(3.0 * t)) }
    private val fD: (Double) -> Double = { t -> 3.0 * cos(3.0 * t) * exp(sin(3.0 * t)) }
    private val fDD: (Double) -> Double = { t ->
        val s = sin(3.0 * t)
        val c = cos(3.0 * t)
        (-9.0 * s + 9.0 * c * c) * exp(s)
    }

    private val systems: List<Pair<String, GeneratingSystem>> = listOf(
        "B" to GeneratingSystem.B,
        "H" to GeneratingSystem.H,
        "T" to GeneratingSystem.T,
    )

    private val families: List<Pair<String, (MinimalSplineBasis) -> FunctionalFamily>> = listOf(
        "theta" to { b -> ProjFunctionals(b) },
        "xi" to { b -> DeBoorFixFunctionals(b, 1) },
        "xitilde" to { b -> DiscreteDeBoorFixFunctionals(b, 1) },
        "mu" to { b -> AveragingFunctionals(b) },
        "lambda" to { b -> ThreePointFunctionals(b) },
    )

    private val grids: List<Pair<String, (Int) -> Grid>> = listOf(
        "uniform" to { n -> Grid.uniform(n, 0.0, 1.0) },
        "quasiUniform" to { n -> Grid.quasiUniform(n, 0.0, 1.0) },
    )

    private val valueWindow = 2.7..3.3
    private val derivWindow = 1.7..2.3

    /** Индексы пар (32→64, 64→128) в списке порядков по [ns]. */
    private val checkedPairs = listOf(2, 3)

    private class Series(val eh: List<Double>) {
        val ord: List<Double> = orders(eh)
    }

    private class Measurement(val value: Series, val interior: Series, val deriv: Series)

    private val results: Map<Triple<String, String, String>, Measurement> by lazy { measureAll() }

    private fun measureAll(): Map<Triple<String, String, String>, Measurement> {
        val out = File("build/reports/convergence-orders.tsv")
        out.parentFile.mkdirs()
        val sb = StringBuilder("sys\tfamily\tgrid\tquantity\tn\tEh\torder\n")
        val map = LinkedHashMap<Triple<String, String, String>, Measurement>()
        for ((sysName, sys) in systems) for ((famName, mk) in families) for ((gridName, mkGrid) in grids) {
            val value = ArrayList<Double>()
            val interior = ArrayList<Double>()
            val deriv = ArrayList<Double>()
            for (n in ns) {
                val grid = mkGrid(n)
                val basis = MinimalSplineBasis(sys, grid)
                val family = mk(basis)
                val gD: (Double) -> Double = if (family.usesDerivative) fD else { _ -> 0.0 }
                val gDD: (Double) -> Double = if (family.usesSecondDerivative) fDD else { _ -> 0.0 }
                val c = family.projectorCoeffs(f, gD, gDD)
                value += errorEh(f, { t -> basis.evalSpline(c, t) }, grid)
                interior += interiorErrorEh(f, { t -> basis.evalSpline(c, t) }, grid)
                deriv += errorEh(fD, { t -> basis.evalSplineDeriv(c, t) }, grid)
            }
            val m = Measurement(Series(value), Series(interior), Series(deriv))
            map[Triple(sysName, famName, gridName)] = m
            for ((q, s) in listOf("value" to m.value, "valueInterior" to m.interior, "deriv" to m.deriv)) {
                for (i in ns.indices) {
                    sb.append("$sysName\t$famName\t$gridName\t$q\t${ns[i]}\t${"%.6e".format(s.eh[i])}\t${"%.4f".format(s.ord[i])}\n")
                }
            }
        }
        out.writeText(sb.toString())
        return map
    }

    private fun fmt(s: Series): String =
        ns.indices.joinToString(", ") { i -> "n=${ns[i]}: E=${"%.3e".format(s.eh[i])} p=${"%.3f".format(s.ord[i])}" }

    private fun assertMonotone(label: String, s: Series) {
        for (i in 0 until s.eh.size - 1) {
            assertTrue(s.eh[i + 1] < s.eh[i], "$label: E_h не убывает при n=${ns[i]}→${ns[i + 1]}; ${fmt(s)}")
        }
    }

    private fun assertOrders(label: String, s: Series, window: ClosedFloatingPointRange<Double>) {
        for (i in checkedPairs) {
            val p = s.ord[i]
            assertTrue(
                p.isFinite() && p in window,
                "$label: порядок на паре n=${ns[i]}→${ns[i + 1]} равен ${"%.4f".format(p)}, окно $window; ${fmt(s)}",
            )
        }
    }

    @TestFactory
    fun convergenceOrders(): List<DynamicTest> {
        val tests = ArrayList<DynamicTest>()
        for ((key, m) in results) {
            val (sysName, famName, gridName) = key
            val label = "$sysName/$famName/$gridName"
            tests += dynamicTest("$label: E_h убывает монотонно") {
                assertMonotone("$label value", m.value)
            }
            if (famName == "xitilde") {
                tests += dynamicTest("$label: порядок 3 по значению во внутренней области") {
                    assertOrders("$label interior", m.interior, valueWindow)
                }
                tests += dynamicTest("$label: порядок 2 по значению с краевым слоем") {
                    assertOrders("$label value", m.value, derivWindow)
                }
            } else {
                tests += dynamicTest("$label: порядок 3 по значению") {
                    assertOrders("$label value", m.value, valueWindow)
                }
            }
            if (famName == "theta" || famName == "xi") {
                tests += dynamicTest("$label: порядок 2 по производной") {
                    assertOrders("$label deriv", m.deriv, derivWindow)
                }
            }
        }
        return tests
    }
}
