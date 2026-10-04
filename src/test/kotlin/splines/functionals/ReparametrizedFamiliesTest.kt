package splines.functionals

import org.junit.jupiter.api.Tag
import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import splines.Reparametrization
import java.io.File
import java.util.Locale
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.pow
import kotlin.test.Test
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Functional families on the reparametrized system G = (1, g, g^2), g(t) = t^beta, over power-graded
 * grids. The de Boor–Fix families must reject G (their coefficients use t-derivatives of the frame,
 * infinite at the left end); the value-based families mu and lambda build on G and reproduce 1, g, g^2.
 * The measured maxima are appended to build/reparametrized-families-measurements.txt.
 */
@Tag("fast")
class ReparametrizedFamiliesTest {
    private val grid: Grid = Grid.power(32, r = 3.0)
    private val betas = listOf(0.5, 1.0 / 3.0)

    private fun basisG(beta: Double): MinimalSplineBasis =
        MinimalSplineBasis(GeneratingSystem.reparametrized(Reparametrization.power(beta)), grid)

    private fun points(g: Grid): List<Double> =
        (0..g.n).map { g.x(it) } + (0 until g.n).flatMap { k -> listOf(0.25, 0.5, 0.75).map { s -> g.x(k) + s * (g.x(k + 1) - g.x(k)) } }

    private fun reproDefect(fam: FunctionalFamily, beta: Double): Double {
        val fs: List<(Double) -> Double> = listOf({ _ -> 1.0 }, { t -> t.pow(beta) }, { t -> t.pow(2.0 * beta) })
        var m = 0.0
        for (f in fs) {
            val c = fam.projectorCoeffs(f)
            for (t in points(fam.grid)) m = max(m, abs(fam.basis.evalSpline(c, t) - f(t)))
        }
        return m
    }

    private fun record(line: String) {
        File("build").mkdirs()
        File("build/reparametrized-families-measurements.txt").appendText(line + "\n")
    }

    @Test
    fun `de Boor-Fix families reject a reparametrized system`() {
        for (beta in betas) {
            val b = basisG(beta)
            for (r in 0..2) {
                val e = assertFailsWith<IllegalArgumentException> { DeBoorFixFunctionals(b, r) }
                assertTrue(e.message!!.contains("de Boor–Fix") && e.message!!.contains(b.sys.name), e.message)
            }
            for (r in 1..2) {
                val e = assertFailsWith<IllegalArgumentException> { DiscreteDeBoorFixFunctionals(b, r) }
                assertTrue(e.message!!.contains("de Boor–Fix") && e.message!!.contains(b.sys.name), e.message)
            }
        }
    }

    @Test
    fun `de Boor-Fix families still build for B on a power grid`() {
        val b = MinimalSplineBasis(GeneratingSystem.B, grid)
        for (r in 0..2) DeBoorFixFunctionals(b, r)
        for (r in 1..2) DiscreteDeBoorFixFunctionals(b, r)
    }

    @Test
    fun `averaging and three-point families on G`() {
        for (beta in betas) {
            val b = basisG(beta)
            val families: List<Pair<String, () -> FunctionalFamily>> = listOf(
                "mu" to { AveragingFunctionals(b) },
                "lambda" to { ThreePointFunctionals(b) },
            )
            for ((name, build) in families) {
                // Both families are value-based and exact on span{1, rho, sigma} = span{1, g, g^2}.
                val defect = reproDefect(build(), beta)
                val line = "$name on ${b.sys.name}, power(32, r=3): max repro defect " +
                    String.format(Locale.ROOT, "%.3e", defect)
                record(line)
                println(line)
                assertTrue(defect <= 1e-12, line)
            }
        }
    }
}
