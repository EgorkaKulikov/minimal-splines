package splines

import org.junit.jupiter.api.DynamicTest
import org.junit.jupiter.api.DynamicTest.dynamicTest
import org.junit.jupiter.api.Tag
import org.junit.jupiter.api.TestFactory
import java.io.File
import kotlin.math.abs
import kotlin.math.max
import kotlin.test.assertTrue

/**
 * An independent oracle for the polynomial case: for φ = (1, t, t²) the minimal splines ω_j equal
 * the classical quadratic B-splines N_{j,2} on the knot vector
 * x_{-2} = x_{-1} = x_0 = a < x_1 < … < x_{n-1} < x_n = x_{n+1} = x_{n+2} = b.
 *
 * N_{j,2} is evaluated by the Cox–de Boor recursion (with the convention 0/0 = 0), the derivative
 * by the formula N'_{j,2} = 2 [N_{j,1}/(x_{j+2} − x_j) − N_{j+1,1}/(x_{j+3} − x_{j+1})]. The oracle
 * implementation uses neither [MinimalSplineBasis] nor [ReferenceSplines].
 *
 * Normalization: both families form a partition of unity, hence the constant of the ratio
 * ω_j / N_{j,2} equals 1; this is checked separately via the ratio inside the support.
 */
@Tag("fast")
class DeBoorOracleTest {
    private val ns = listOf(4, 8, 16)
    private val relTol = 1e-13

    private val grids: List<Pair<String, (Int) -> Grid>> = listOf(
        "uniform" to { n -> Grid.uniform(n, 0.0, 1.0) },
        "quasiUniform" to { n -> Grid.quasiUniform(n, 0.0, 1.0) },
        "geometric" to { n -> Grid.geometric(n, 0.0, 1.0) },
        "graded" to { n -> Grid.graded(n, 0.0, 1.0) },
    )

    /** N_{j,k}(t) by Cox–de Boor; at t = b the last non-degenerate interval is treated as closed. */
    private fun bspline(grid: Grid, j: Int, k: Int, t: Double): Double {
        if (k == 0) {
            val xj = grid.x(j)
            val xj1 = grid.x(j + 1)
            if (t == grid.b) return if (xj < xj1 && xj1 == grid.b) 1.0 else 0.0
            return if (xj <= t && t < xj1) 1.0 else 0.0
        }
        val left = grid.x(j + k) - grid.x(j)
        val right = grid.x(j + k + 1) - grid.x(j + 1)
        val w1 = if (left > 0.0) (t - grid.x(j)) / left else 0.0
        val w2 = if (right > 0.0) (grid.x(j + k + 1) - t) / right else 0.0
        return w1 * bspline(grid, j, k - 1, t) + w2 * bspline(grid, j + 1, k - 1, t)
    }

    private fun bsplineDeriv(grid: Grid, j: Int, t: Double): Double {
        val d1 = grid.x(j + 2) - grid.x(j)
        val d2 = grid.x(j + 3) - grid.x(j + 1)
        val a = if (d1 > 0.0) bspline(grid, j, 1, t) / d1 else 0.0
        val b = if (d2 > 0.0) bspline(grid, j + 1, 1, t) / d2 else 0.0
        return 2.0 * (a - b)
    }

    private fun controlPoints(grid: Grid): DoubleArray {
        val m = 20 * grid.n
        val pts = DoubleArray(m + 1) { i -> grid.a + (grid.b - grid.a) * i / m }
        return (pts.toList() + (0..grid.n).map { grid.x(it) }).distinct().sorted().toDoubleArray()
    }

    private class Record(
        val maxOmega: Double, val maxDiffOmega: Double, val whereOmega: String,
        val maxDeriv: Double, val maxDiffDeriv: Double, val whereDeriv: String,
        val maxRatioDev: Double,
    )

    /** All measurements; as a side effect `build/reports/deboor-oracle.tsv` is written. */
    private val results: Map<Pair<String, Int>, Record> by lazy { measureAll() }

    private fun measureAll(): Map<Pair<String, Int>, Record> {
        val out = File("build/reports/deboor-oracle.tsv")
        out.parentFile.mkdirs()
        val sb = StringBuilder("grid\tn\tmaxN\tmaxDiffOmega\tmaxDN\tmaxDiffOmegaDeriv\tmaxRatioDev\n")
        val map = LinkedHashMap<Pair<String, Int>, Record>()
        for ((gridName, mkGrid) in grids) for (n in ns) {
            val grid = mkGrid(n)
            val basis = MinimalSplineBasis(GeneratingSystem.B, grid)
            val pts = controlPoints(grid)
            var maxOmega = 0.0
            var maxDiffOmega = 0.0
            var whereOmega = ""
            var maxDeriv = 0.0
            var maxDiffDeriv = 0.0
            var whereDeriv = ""
            for (j in -2..n - 1) for (t in pts) {
                val ref = bspline(grid, j, 2, t)
                maxOmega = max(maxOmega, abs(ref))
                val d = abs(basis.omega(j, t) - ref)
                if (d > maxDiffOmega) {
                    maxDiffOmega = d
                    whereOmega = "j=$j t=$t omega=${basis.omega(j, t)} N=$ref"
                }
                val refD = bsplineDeriv(grid, j, t)
                maxDeriv = max(maxDeriv, abs(refD))
                val dD = abs(basis.omegaDeriv(j, t) - refD)
                if (dD > maxDiffDeriv) {
                    maxDiffDeriv = dD
                    whereDeriv = "j=$j t=$t omega'=${basis.omegaDeriv(j, t)} N'=$refD"
                }
            }
            var maxRatioDev = 0.0
            for (j in -2..n - 1) {
                val lo = grid.x(j)
                val hi = grid.x(j + 3)
                for (q in 1..7) {
                    val t = lo + (hi - lo) * q / 8.0
                    val ref = bspline(grid, j, 2, t)
                    if (ref < 1e-3) continue
                    maxRatioDev = max(maxRatioDev, abs(basis.omega(j, t) / ref - 1.0))
                }
            }
            map[gridName to n] = Record(maxOmega, maxDiffOmega, whereOmega, maxDeriv, maxDiffDeriv, whereDeriv, maxRatioDev)
            sb.append("$gridName\t$n\t${"%.6e".format(maxOmega)}\t${"%.3e".format(maxDiffOmega)}\t${"%.6e".format(maxDeriv)}\t${"%.3e".format(maxDiffDeriv)}\t${"%.3e".format(maxRatioDev)}\n")
        }
        out.writeText(sb.toString())
        return map
    }

    @TestFactory
    fun omegaMatchesCoxDeBoor(): List<DynamicTest> {
        val tests = ArrayList<DynamicTest>()
        for ((key, r) in results) {
            val (gridName, n) = key
            val label = "B/$gridName/n=$n"
            tests += dynamicTest("$label: omega_j = N_{j,2}") {
                assertTrue(
                    r.maxDiffOmega <= relTol * r.maxOmega,
                    "$label: max|omega_j − N_{j,2}| = ${r.maxDiffOmega} (max|N| = ${r.maxOmega}) at ${r.whereOmega}",
                )
            }
            tests += dynamicTest("$label: omega_j' = N'_{j,2}") {
                assertTrue(
                    r.maxDiffDeriv <= relTol * r.maxDeriv,
                    "$label: max|omega_j' − N'_{j,2}| = ${r.maxDiffDeriv} (max|N'| = ${r.maxDeriv}) at ${r.whereDeriv}",
                )
            }
            tests += dynamicTest("$label: normalization constant omega_j / N_{j,2} equals 1") {
                assertTrue(r.maxRatioDev <= 1e-12, "$label: max|omega_j / N_{j,2} − 1| = ${r.maxRatioDev}")
            }
        }
        return tests
    }
}
