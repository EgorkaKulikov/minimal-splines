package splines.functionals

import numerics.NumericsContext
import org.junit.jupiter.api.Tag
import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import splines.Reparametrization
import java.io.File
import java.util.Locale
import kotlin.math.abs
import kotlin.math.max
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * [ThetaSampling] for [ProjFunctionals] on the systems B and G = reparametrized(power(beta)):
 * biorthogonality and reproduction for both samplings, the equivalence of theta^tau on G with
 * theta for B over the tau-grid, the constants C_chi, and the compatibility of the default
 * constructor. The measured maxima and the C_chi table are appended to
 * build/theta-sampling-measurements.txt.
 */
@Tag("fast")
class ProjFunctionalsSamplingTest {
    private val betas = listOf(0.5, 2.0 / 3.0)
    private val ns = listOf(8, 32, 128)
    private val grids: List<Pair<String, (Int) -> Grid>> = listOf(
        "uniform" to { n -> Grid.uniform(n) },
        "power2" to { n -> Grid.power(n, r = 2.0) },
        "power3" to { n -> Grid.power(n, r = 3.0) },
    )
    private val fractions = doubleArrayOf(0.0, 0.07, 0.31, 0.5, 0.77, 0.98)
    private val tol = 1e-12
    private val out = File("build/theta-sampling-measurements.txt")

    private data class Case(val beta: Double, val gridName: String, val grid: Grid) {
        val rep: Reparametrization = Reparametrization.power(beta)
        val tauSampling: ThetaSampling = ThetaSampling.reparametrized(rep)
        val basisG: MinimalSplineBasis = MinimalSplineBasis(GeneratingSystem.reparametrized(rep), grid)
        val basisB: MinimalSplineBasis = MinimalSplineBasis(GeneratingSystem.B, grid)
        override fun toString(): String = "beta=${fmt(beta)}, $gridName, n=${grid.n}"
    }

    private val cases: List<Case> by lazy {
        betas.flatMap { beta -> grids.flatMap { (name, make) -> ns.map { n -> Case(beta, name, make(n)) } } }
    }

    private fun points(grid: Grid): List<Double> =
        (0 until grid.n).flatMap { k ->
            val l = grid.x(k); val r = grid.x(k + 1)
            fractions.map { l + it * (r - l) }
        } + grid.b

    private fun bioDefect(theta: ProjFunctionals): Double {
        var m = 0.0
        for (i in -2..theta.n - 1) for (j in -2..theta.n - 1) {
            val v = theta.chi(i).apply { t -> theta.basis.omega(j, t) }
            m = max(m, abs(v - if (i == j) 1.0 else 0.0))
        }
        return m
    }

    private fun reproDefect(theta: ProjFunctionals, fs: List<(Double) -> Double>): Double {
        var m = 0.0
        for (f in fs) {
            val coeffs = theta.projectorCoeffs(f)
            for (t in points(theta.grid)) m = max(m, abs(theta.basis.evalSpline(coeffs, t) - f(t)))
        }
        return m
    }

    private fun value(theta: ProjFunctionals, j: Int): ValueFunctional = theta.chi(j) as ValueFunctional

    @Test
    fun `biorthogonality and reproduction for both samplings on B and G`() {
        var bio = 0.0; var repro = 0.0
        for (c in cases) {
            val g = c.rep.g
            val fsG: List<(Double) -> Double> = listOf({ 1.0 }, g, { t -> g(t) * g(t) })
            val fsB: List<(Double) -> Double> = listOf({ 1.0 }, { t -> t }, { t -> t * t })
            for (sampling in listOf(ThetaSampling.ARITHMETIC, c.tauSampling)) {
                val thG = ProjFunctionals(c.basisG, NumericsContext.default(), sampling)
                val thB = ProjFunctionals(c.basisB, NumericsContext.default(), sampling)
                val b = max(bioDefect(thG), bioDefect(thB))
                val r = max(reproDefect(thG, fsG), reproDefect(thB, fsB))
                assertTrue(b <= tol, "biorthogonality defect $b for $c, ${thG.name}")
                assertTrue(r <= tol, "reproduction defect $r for $c, ${thG.name}")
                bio = max(bio, b); repro = max(repro, r)
            }
        }
        out.parentFile.mkdirs()
        out.appendText("ThetaSampling: max biorthogonality defect = $bio, max reproduction defect = $repro\n")
    }

    @Test
    fun `theta-tau on G coincides with theta for B over the tau-grid`() {
        var coeffDev = 0.0; var nodeDev = 0.0; var lambdaDev = 0.0
        val table = StringBuilder("ThetaSampling: C_chi table (B on X, G-t, G-tau, B on tau-grid)\n")
        for (c in cases) {
            val n = c.grid.n
            val tauGrid = Grid(n, DoubleArray(n + 1) { c.rep.g(c.grid.x(it)) })
            val thTau = ProjFunctionals(c.basisG, NumericsContext.default(), c.tauSampling)
            val thBTau = ProjFunctionals(MinimalSplineBasis(GeneratingSystem.B, tauGrid))
            for (j in -2..n - 1) {
                val a = value(thTau, j); val b = value(thBTau, j)
                assertEquals(b.nodes.size, a.nodes.size, "node count at j=$j for $c")
                for (k in a.nodes.indices) {
                    val dc = abs(a.coeffs[k] - b.coeffs[k]) / max(1.0, abs(b.coeffs[k]))
                    val back = c.rep.gInverse(b.nodes[k])
                    val dn = abs(a.nodes[k] - back) / max(abs(a.nodes[k]), Double.MIN_VALUE)
                    assertTrue(dc <= tol, "coefficient deviation $dc at j=$j, k=$k for $c")
                    assertTrue(dn <= tol, "node deviation $dn at j=$j, k=$k for $c")
                    coeffDev = max(coeffDev, dc); nodeDev = max(nodeDev, dn)
                }
            }
            val lamTau = thTau.cChi(); val lamBTau = thBTau.cChi()
            val dl = abs(lamTau - lamBTau) / lamBTau
            assertTrue(dl <= tol, "C_chi deviation $dl for $c")
            lambdaDev = max(lambdaDev, dl)
            val lamB = ProjFunctionals(c.basisB).cChi()
            val lamGt = ProjFunctionals(c.basisG).cChi()
            table.append("  $c: B=${fmt(lamB)} G-t=${fmt(lamGt)} G-tau=${fmt(lamTau)} B(tau)=${fmt(lamBTau)}\n")
        }
        out.parentFile.mkdirs()
        out.appendText(table.toString())
        out.appendText(
            "ThetaSampling: G-tau vs B(tau-grid): max coefficient rel dev = $coeffDev, " +
                "max node rel dev = $nodeDev, max C_chi rel dev = $lambdaDev\n",
        )
    }

    @Test
    fun `default constructor keeps the name theta and the arithmetic coefficients bitwise`() {
        for (c in cases) {
            val default = ProjFunctionals(c.basisG)
            val withCtx = ProjFunctionals(c.basisG, NumericsContext.default())
            val explicit = ProjFunctionals(c.basisG, NumericsContext.default(), ThetaSampling.ARITHMETIC)
            val tau = ProjFunctionals(c.basisG, NumericsContext.default(), c.tauSampling)
            assertEquals("theta", default.name)
            assertEquals("theta", withCtx.name)
            assertEquals("theta", explicit.name)
            assertEquals("theta-tau", tau.name)
            assertTrue(default.sampling === ThetaSampling.ARITHMETIC)
            for (j in -2..c.grid.n - 1) {
                assertContentEquals(value(explicit, j).nodes, value(default, j).nodes, "nodes j=$j, $c")
                assertContentEquals(value(explicit, j).coeffs, value(default, j).coeffs, "coeffs j=$j, $c")
                assertContentEquals(value(explicit, j).coeffs, value(withCtx, j).coeffs, "coeffs j=$j, $c")
            }
        }
    }

    private companion object {
        fun fmt(v: Double): String = String.format(Locale.ROOT, "%.6f", v)
    }
}
