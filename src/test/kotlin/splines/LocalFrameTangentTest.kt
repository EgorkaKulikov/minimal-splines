package splines

import org.junit.jupiter.api.Tag
import kotlin.math.abs
import kotlin.math.max
import kotlin.test.Test
import kotlin.test.assertSame
import kotlin.test.assertTrue

/**
 * Contract of [LocalFrame.tangent]: it equals psiD for the built-in systems, and the basis does not
 * depend on the finite non-vanishing factor between tangent and psiD, because the approximation-relation
 * vectors a_j are invariant to that factor.
 */
@Tag("fast")
class LocalFrameTangentTest {
    private val systems = listOf(GeneratingSystem.B, GeneratingSystem.H, GeneratingSystem.T)

    private val factors: List<Pair<String, (Double) -> Double>> = listOf(
        "1 + t^2" to { t -> 1.0 + t * t },
        "1/(1 + t)" to { t -> 1.0 / (1.0 + t) },
    )

    private fun grids(): List<Grid> = listOf(8, 16).flatMap { n ->
        listOf(Grid.geometric(n, R = 4.0), Grid.quasiUniform(n))
    }

    /** The same system as [sys] whose local frames carry tangent = lambda(t) psiD(t). */
    private fun scaledTangentSystem(sys: GeneratingSystem, lambda: (Double) -> Double): GeneratingSystem =
        GeneratingSystem(
            name = sys.name + "~",
            rho = sys.rho, sigma = sys.sigma,
            rhoD = sys.rhoD, sigmaD = sys.sigmaD,
            rhoDD = sys.rhoDD, sigmaDD = sys.sigmaDD,
            localFrameFactory = { c, h ->
                val f = sys.localFrame(c, h)
                LocalFrame(f.psi, f.psiD, f.psiDD, f.det) { t ->
                    val d = f.psiD(t); val l = lambda(t)
                    doubleArrayOf(l * d[0], l * d[1], l * d[2])
                }
            },
        )

    @Test fun builtInFramesUsePsiDAsTangent() {
        for (sys in systems) {
            for ((c, h) in listOf(0.0 to 0.1, 0.3 to 0.05, 2.0 to 1.5)) {
                val f = sys.localFrame(c, h)
                assertSame(f.psiD, f.tangent, "tangent !== psiD for ${sys.name} at c=$c, h=$h")
            }
        }
    }

    @Test fun globalFrameUsesPhiDAsTangent() {
        val custom = GeneratingSystem(
            name = "B-global",
            rho = { t -> t }, sigma = { t -> t * t },
            rhoD = { 1.0 }, sigmaD = { t -> 2.0 * t },
            rhoDD = { 0.0 }, sigmaDD = { 2.0 },
        )
        val f = custom.localFrame(0.2, 0.1)
        assertSame(f.psiD, f.tangent)
    }

    @Test fun fourArgumentConstructorSetsTangentToPsiD() {
        val psi: (Double) -> DoubleArray = { t -> doubleArrayOf(1.0, t, t * t) }
        val psiD: (Double) -> DoubleArray = { t -> doubleArrayOf(0.0, 1.0, 2.0 * t) }
        val psiDD: (Double) -> DoubleArray = { doubleArrayOf(0.0, 0.0, 2.0) }
        val f = LocalFrame(psi, psiD, psiDD, 1.0)
        assertSame(psiD, f.tangent)
        assertSame(psiD, f.psiD)
    }

    @Test fun computeAIsInvariantToTangentFactor() {
        var worst = 0.0
        for (sys in systems) for (grid in grids()) {
            val basis = MinimalSplineBasis(sys, grid)
            for ((label, lambda) in factors) {
                for (j in -2..grid.n - 1) {
                    for (k in max(0, j)..minOf(grid.n - 1, j + 2)) {
                        val f = basis.frame(k)
                        val ref = basis.computeA(j, f.psi, f.psiD)
                        val scaled = basis.computeA(j, f.psi) { t ->
                            val d = f.psiD(t); val l = lambda(t)
                            doubleArrayOf(l * d[0], l * d[1], l * d[2])
                        }
                        val err = relDiff(ref, scaled)
                        worst = max(worst, err)
                        assertTrue(err <= 1e-14, "${sys.name}, n=${grid.n}, lambda=$label, j=$j, k=$k: rel diff $err")
                    }
                }
            }
        }
        println("LocalFrameTangentTest.computeAIsInvariantToTangentFactor: max rel diff = $worst")
    }

    @Test fun basisIsInvariantToTangentFactor() {
        var worst = 0.0
        for (sys in systems) for (grid in grids()) {
            val ref = MinimalSplineBasis(sys, grid)
            for ((label, lambda) in factors) {
                val scaled = MinimalSplineBasis(scaledTangentSystem(sys, lambda), grid)
                val m = 40 * grid.n
                for (i in 0..m) {
                    val t = grid.a + (grid.b - grid.a) * i / m
                    for (j in -2..grid.n - 1) {
                        val r = ref.omega(j, t)
                        val s = scaled.omega(j, t)
                        val err = abs(s - r) / max(1.0, abs(r))
                        worst = max(worst, err)
                        assertTrue(err <= 1e-13, "${sys.name}, n=${grid.n}, lambda=$label, j=$j, t=$t: rel diff $err")
                    }
                }
            }
        }
        println("LocalFrameTangentTest.basisIsInvariantToTangentFactor: max rel diff = $worst")
    }

    private fun relDiff(a: DoubleArray, b: DoubleArray): Double {
        var num = 0.0; var den = 0.0
        for (i in a.indices) { num = max(num, abs(a[i] - b[i])); den = max(den, abs(a[i])) }
        return num / max(den, Double.MIN_VALUE)
    }
}
