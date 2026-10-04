package splines

import numerics.Conditioning
import numerics.backend.Backends
import org.junit.jupiter.api.Tag
import splines.functionals.ProjFunctionals
import splines.functionals.apply
import java.io.File
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.pow
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertSame
import kotlin.test.assertTrue

/**
 * The reparametrized system phi = (1, g, g^2), g(t) = t^beta: its minimal splines are
 * omega_j = B_j ∘ g, where B_j are the quadratic B-splines on the image knots g(x_j) (triple at the
 * ends, as the grid knots). The oracle evaluates B_j by the Cox–de Boor recursion and uses neither
 * [MinimalSplineBasis] nor [ReferenceSplines]. The measured maxima are appended to
 * build/reparametrized-measurements.txt.
 */
@Tag("fast")
class ReparametrizedSystemTest {
    private val betas = listOf(0.5, 2.0 / 3.0, 1.0 / 3.0)
    private val ns = listOf(8, 16, 64)
    private val grids: List<Pair<String, (Int) -> Grid>> = listOf(
        "uniform" to { n -> Grid.uniform(n) },
        "geometric" to { n -> Grid.geometric(n) },
    )
    private val fractions = doubleArrayOf(0.0, 0.03, 0.2, 0.45, 0.7, 0.9, 0.99)
    private val backend = Backends.default()

    private data class Case(val beta: Double, val gridName: String, val grid: Grid) {
        val rep: Reparametrization = Reparametrization.power(beta)
        val basis: MinimalSplineBasis = MinimalSplineBasis(GeneratingSystem.reparametrized(rep), grid)
        override fun toString(): String = "beta=$beta, $gridName, n=${grid.n}"
    }

    private val cases: List<Case> by lazy {
        betas.flatMap { beta -> grids.flatMap { (name, make) -> ns.map { n -> Case(beta, name, make(n)) } } }
    }

    private fun points(grid: Grid): DoubleArray {
        val inner = (0 until grid.n).flatMap { k ->
            val l = grid.x(k); val r = grid.x(k + 1)
            fractions.map { f -> l + f * (r - l) }
        }
        return (listOf(1e-12, 1e-8) + inner + listOf(grid.b)).toDoubleArray()
    }

    /** Image knots y_{j} = g(x_j), j = -2..n+2, stored with offset 2. */
    private fun imageKnots(grid: Grid, g: (Double) -> Double): DoubleArray =
        DoubleArray(grid.n + 5) { g(grid.x(it - 2)) }

    /** N_{j,k}(v) on the knot vector y by Cox–de Boor (0/0 = 0); at v = y_last the last interval is closed. */
    private fun bspline(y: DoubleArray, j: Int, k: Int, v: Double): Double {
        val i = j + 2
        if (k == 0) {
            val l = y[i]; val r = y[i + 1]
            if (v == y.last()) return if (l < r && r == y.last()) 1.0 else 0.0
            return if (l <= v && v < r) 1.0 else 0.0
        }
        val left = y[i + k] - y[i]
        val right = y[i + k + 1] - y[i + 1]
        val w1 = if (left > 0.0) (v - y[i]) / left else 0.0
        val w2 = if (right > 0.0) (y[i + k + 1] - v) / right else 0.0
        return w1 * bspline(y, j, k - 1, v) + w2 * bspline(y, j + 1, k - 1, v)
    }

    /** N'_{j,2}(v) = 2 [N_{j,1}/(y_{j+2} - y_j) - N_{j+1,1}/(y_{j+3} - y_{j+1})]. */
    private fun bsplineDeriv(y: DoubleArray, j: Int, v: Double): Double {
        val i = j + 2
        val d1 = y[i + 2] - y[i]
        val d2 = y[i + 3] - y[i + 1]
        val a = if (d1 > 0.0) bspline(y, j, 1, v) / d1 else 0.0
        val b = if (d2 > 0.0) bspline(y, j + 1, 1, v) / d2 else 0.0
        return 2.0 * (a - b)
    }

    private fun record(line: String) {
        File("build").mkdirs()
        File("build/reparametrized-measurements.txt").appendText(line + "\n")
    }

    @Test
    fun `omega_j equals B_j composed with g`() {
        var worst = 0.0
        for (c in cases) {
            val y = imageKnots(c.grid, c.rep.g)
            var mx = 0.0
            for (t in points(c.grid)) {
                val v = c.rep.g(t)
                for (j in -2..c.grid.n - 1) mx = max(mx, abs(c.basis.omega(j, t) - bspline(y, j, 2, v)))
            }
            assertTrue(mx <= 1e-12, "omega_j != B_j ∘ g for $c: max error $mx")
            worst = max(worst, mx)
        }
        record("B∘g agreement max abs error: $worst")
    }

    @Test
    fun `beta equal to one reproduces the polynomial basis`() {
        var worst = 0.0
        for ((_, make) in grids) for (n in ns) {
            val grid = make(n)
            val rep = MinimalSplineBasis(GeneratingSystem.reparametrized(Reparametrization.power(1.0)), grid)
            val ref = MinimalSplineBasis(GeneratingSystem.B, grid)
            for (t in points(grid)) for (j in -2..n - 1) {
                worst = max(worst, abs(rep.omega(j, t) - ref.omega(j, t)))
                if (t > 0.0 && t < grid.b) {
                    worst = max(worst, abs(rep.omegaDeriv(j, t) - ref.omegaDeriv(j, t)) / n)
                    worst = max(worst, abs(rep.omegaDeriv2(j, t) - ref.omegaDeriv2(j, t)) / (n * n))
                }
            }
        }
        assertTrue(worst <= 1e-13, "beta = 1 differs from B: $worst")
        record("beta=1 vs B max abs error (derivatives scaled by h, h^2): $worst")
    }

    @Test
    fun `partition of unity, biorthogonality and reproduction of 1, g, g^2`() {
        var pu = 0.0; var bio = 0.0; var repro = 0.0
        for (c in cases) {
            val n = c.grid.n
            val theta = ProjFunctionals(c.basis)
            for (t in points(c.grid)) pu = max(pu, abs((-2..n - 1).sumOf { c.basis.omega(it, t) } - 1.0))
            for (i in -2..n - 1) for (j in -2..n - 1) {
                val v = theta.chi(i).apply { t -> c.basis.omega(j, t) }
                bio = max(bio, abs(v - if (i == j) 1.0 else 0.0))
            }
            val g = c.rep.g
            val fs: List<(Double) -> Double> = listOf({ 1.0 }, g, { t -> g(t) * g(t) })
            for (f in fs) {
                val coeffs = theta.projectorCoeffs(f)
                for (t in points(c.grid)) repro = max(repro, abs(c.basis.evalSpline(coeffs, t) - f(t)))
            }
        }
        assertTrue(pu <= 1e-12, "partition of unity: $pu")
        assertTrue(bio <= 1e-12, "biorthogonality: $bio")
        assertTrue(repro <= 1e-12, "reproduction of 1, g, g^2: $repro")
        record("partition of unity: $pu; biorthogonality: $bio; reproduction: $repro")
    }

    @Test
    fun `local approximation matrices are well conditioned`() {
        var worst = 0.0
        for (c in cases) for (k in 0 until c.grid.n) {
            val cond = Conditioning.conditionEstimate(c.basis.localApproximationMatrix(k), backend = backend).condInf
            assertTrue(cond <= MinimalSplineBasis.MAX_CONDITION, "cond $cond on interval $k for $c")
            worst = max(worst, cond)
        }
        record("max condition estimate of T_k M_k: $worst")
    }

    @Test
    fun `derivatives follow the chain rule`() {
        val c = Case(0.5, "uniform", Grid.uniform(8))
        val y = imageKnots(c.grid, c.rep.g)
        for (k in 0 until c.grid.n) {
            val h = c.grid.x(k + 1) - c.grid.x(k)
            val t = c.grid.x(k) + 0.5 * h
            for (j in k - 2..k) {
                val d1 = bsplineDeriv(y, j, c.rep.g(t)) * c.rep.gD(t)
                assertEquals(d1, c.basis.omegaDeriv(j, t), 1e-10 * max(1.0, abs(d1)), "omega' j=$j t=$t")
                val s = 1e-5 * h
                val fd = (c.basis.omegaDeriv(j, t + s) - c.basis.omegaDeriv(j, t - s)) / (2.0 * s)
                assertEquals(fd, c.basis.omegaDeriv2(j, t), 1e-5 * max(1.0, abs(fd)), "omega'' j=$j t=$t")
            }
        }
    }

    @Test
    fun `system properties, wronskian and inverse map`() {
        val rep = Reparametrization.power(0.5)
        val sys = GeneratingSystem.reparametrized(rep)
        assertEquals("G[pow(0.5)]", sys.name)
        assertSame(rep, sys.reparametrization)
        assertNull(GeneratingSystem.B.reparametrization)
        assertNull(GeneratingSystem(GeneratingSystem.H.name, { it }, { it }, { 1.0 }, { 1.0 }, { 0.0 }, { 0.0 }).reparametrization)
        for (t in doubleArrayOf(0.01, 0.3, 0.9)) {
            val w = 2.0 * rep.gD(t).pow(3)
            assertEquals(w, sys.wronskian(t), 1e-12 * w, "wronskian at t=$t")
            assertEquals(t, rep.gInverse(rep.g(t)), 1e-15, "inverse at t=$t")
        }
        assertTrue(rep.difference(0.0, 0.0) == 0.0)
        val shifted = Reparametrization.power(1.0 / 3.0, a = -1.0)
        assertEquals(shifted.g(0.5) - shifted.g(0.25), shifted.difference(0.5, 0.25), 1e-15)
    }

    @Test
    fun `default difference gives the same basis`() {
        val linear = Reparametrization("id", { it }, { 1.0 }, { 0.0 }, { it })
        val grid = Grid.geometric(16)
        val a = MinimalSplineBasis(GeneratingSystem.reparametrized(linear), grid)
        val b = MinimalSplineBasis(GeneratingSystem.B, grid)
        for (t in points(grid)) for (j in -2..15) assertEquals(b.omega(j, t), a.omega(j, t), 1e-13)
    }

    @Test
    fun `power validates its arguments`() {
        for (beta in doubleArrayOf(0.0, -0.5, 1.5, Double.NaN)) {
            assertFailsWith<IllegalArgumentException>("beta=$beta") { Reparametrization.power(beta) }
        }
        assertFailsWith<IllegalArgumentException> { Reparametrization.power(0.5, a = Double.NaN) }
    }
}
