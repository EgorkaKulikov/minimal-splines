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
 * Convergence orders of quasi-interpolants on quadratic minimal splines under grid refinement:
 * n = 8, 16, 32, 64, 128 on [0, 1], with the test function f(t) = exp(sin 3t), which lies in
 * span φ of none of the systems B, H, T.
 *
 * Checked:
 *  - monotone decrease of E_h;
 *  - the order by value on the pairs 32→64 and 64→128 within the window [2.7, 3.3] for θ, ξ, μ, λ
 *    (theoretical order 3);
 *  - for ξ̃ — order 3 in the interior region (excluding the boundary layer of width 3h) and order 2
 *    in the full norm: the one-sided difference at multiple boundary nodes produces a boundary
 *    defect O(h²);
 *  - the order by derivative for θ and ξ within the window [1.7, 2.3] (theoretical order 2).
 *
 * All measurements are written to `build/reports/convergence-orders.tsv`
 * (columns `sys family grid quantity n Eh order`).
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

    /** Indices of the pairs (32→64, 64→128) in the list of orders over [ns]. */
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
            assertTrue(s.eh[i + 1] < s.eh[i], "$label: E_h does not decrease at n=${ns[i]}→${ns[i + 1]}; ${fmt(s)}")
        }
    }

    private fun assertOrders(label: String, s: Series, window: ClosedFloatingPointRange<Double>) {
        for (i in checkedPairs) {
            val p = s.ord[i]
            assertTrue(
                p.isFinite() && p in window,
                "$label: the order on the pair n=${ns[i]}→${ns[i + 1]} is ${"%.4f".format(p)}, window $window; ${fmt(s)}",
            )
        }
    }

    @TestFactory
    fun convergenceOrders(): List<DynamicTest> {
        val tests = ArrayList<DynamicTest>()
        for ((key, m) in results) {
            val (sysName, famName, gridName) = key
            val label = "$sysName/$famName/$gridName"
            tests += dynamicTest("$label: E_h decreases monotonically") {
                assertMonotone("$label value", m.value)
            }
            if (famName == "xitilde") {
                tests += dynamicTest("$label: order 3 by value in the interior region") {
                    assertOrders("$label interior", m.interior, valueWindow)
                }
                tests += dynamicTest("$label: order 2 by value including the boundary layer") {
                    assertOrders("$label value", m.value, derivWindow)
                }
            } else {
                tests += dynamicTest("$label: order 3 by value") {
                    assertOrders("$label value", m.value, valueWindow)
                }
            }
            if (famName == "theta" || famName == "xi") {
                tests += dynamicTest("$label: order 2 by derivative") {
                    assertOrders("$label deriv", m.deriv, derivWindow)
                }
            }
        }
        return tests
    }
}
