package splines

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
import kotlin.test.assertTrue

/** Power-graded grids [Grid.power] and [Grid.symmetricPower]: nodes, neighbour ratios, guards, bases on them. */
@Tag("fast")
class GridPowerTest {
    private val exponents = listOf(1.5, 2.0, 3.0, 4.0, 6.0)
    private val sizes = listOf(8, 64, 256)

    private fun record(line: String) {
        File("build").mkdirs()
        File("build/grid-power-measurements.txt").appendText(line + "\n")
    }

    private fun nodes(g: Grid): DoubleArray = DoubleArray(g.n + 1) { g.x(it) }

    private fun isStrictlyIncreasing(x: DoubleArray): Boolean = (0 until x.size - 1).all { x[it] < x[it + 1] }

    /** Max over nodes j in [js] of max(h_j/h_{j-1}, h_{j-1}/h_j) and the index where it is attained. */
    private fun maxNeighbourRatio(x: DoubleArray, js: IntRange = 1 until x.size - 1): Pair<Double, Int> {
        var best = 0.0
        var arg = -1
        for (j in js) {
            val hp = x[j] - x[j - 1]
            val hj = x[j + 1] - x[j]
            val q = max(hj / hp, hp / hj)
            if (q > best) { best = q; arg = j }
        }
        return best to arg
    }

    @Test
    fun `power grid has exact ends, increases strictly and reduces to uniform for r = 1`() {
        for ((a, b) in listOf(0.0 to 1.0, 0.0 to 3.0, -1.0 to 1.0, 2.0 to 3.0, -3.5 to 0.5)) for (n in listOf(1, 2, 7, 64)) {
            for (r in exponents) {
                val x = nodes(Grid.power(n, a, b, r))
                assertEquals(a, x[0]); assertEquals(b, x[n])
                assertTrue(isStrictlyIncreasing(x), "power n=$n r=$r on [$a, $b]")
            }
            val p = Grid.power(n, a, b, 1.0)
            val u = Grid.uniform(n, a, b)
            for (j in -2..n + 2) {
                assertTrue(abs(p.x(j) - u.x(j)) <= Math.ulp(u.x(j)), "r=1 vs uniform: j=$j n=$n on [$a, $b]")
            }
        }
    }

    @Test
    fun `power grid neighbour ratio is 2^r - 1 attained at j = 1 for every n`() {
        for (r in exponents) for (n in sizes) {
            val expected = 2.0.pow(r) - 1.0
            if (n.toDouble().pow(-r) <= 4.0 * Grid.BREAKPOINT_INCLUSION_EPS_UNIT) {
                assertFailsWith<IllegalArgumentException> { Grid.power(n, r = r) }
                record("power r=$r n=$n: rejected by the first-step guard (x_1 = ${n.toDouble().pow(-r)})")
                continue
            }
            val (q, arg) = maxNeighbourRatio(nodes(Grid.power(n, r = r)))
            val rel = abs(q - expected) / expected
            record("power r=$r n=$n: max ratio=$q at j=$arg, 2^r-1=$expected, rel=$rel")
            assertTrue(rel <= 1e-12, "power r=$r n=$n: max ratio $q vs $expected (rel $rel)")
            assertEquals(1, arg, "power r=$r n=$n: argmax")
        }
    }

    @Test
    fun `symmetric power grid is mirror symmetric with neighbour ratio 2^r - 1`() {
        for ((a, b) in listOf(0.0 to 1.0, -1.0 to 1.0, 2.0 to 3.0)) for (n in listOf(2, 8, 64)) for (r in exponents) {
            val x = nodes(Grid.symmetricPower(n, a, b, r))
            val scale = max(abs(a), abs(b))
            assertEquals(a, x[0]); assertEquals(b, x[n]); assertEquals((a + b) / 2, x[n / 2])
            assertTrue(isStrictlyIncreasing(x), "symmetricPower n=$n r=$r on [$a, $b]")
            for (j in 0..n) {
                assertTrue(abs(x[j] + x[n - j] - (a + b)) <= 2 * Math.ulp(scale), "mirror: j=$j n=$n r=$r on [$a, $b]")
            }
        }
        for (r in exponents) for (n in sizes) {
            val expected = 2.0.pow(r) - 1.0
            val x = nodes(Grid.symmetricPower(n, r = r))
            val m = n / 2
            // Left half: computed from the power law directly, so the ratio is exact to rounding of pow.
            val (qLeft, argLeft) = maxNeighbourRatio(x, 1..m)
            // Whole grid: the steps near b are differences of nodes of size 1, so their ratios carry a
            // rounding error of order ulp(1)/h_0, documented in the KDoc of symmetricPower.
            val (q, arg) = maxNeighbourRatio(x)
            val h0 = x[1] - x[0]
            val relLeft = abs(qLeft - expected) / expected
            val rel = abs(q - expected) / expected
            record("symmetricPower r=$r n=$n: left max ratio=$qLeft at j=$argLeft (rel $relLeft); " +
                "full max ratio=$q at j=$arg (rel $rel), h_0=$h0")
            assertTrue(relLeft <= 1e-12, "symmetricPower r=$r n=$n: left max ratio $qLeft vs $expected")
            assertEquals(1, argLeft, "symmetricPower r=$r n=$n: left argmax")
            assertTrue(rel <= 1e-12 + 8 * Math.ulp(1.0) / h0, "symmetricPower r=$r n=$n: max ratio $q vs $expected")
            assertTrue(arg == 1 || arg == n - 1, "symmetricPower r=$r n=$n: argmax $arg")
        }
    }

    @Test
    fun `first-step guard rejects steps not exceeding four breakpoint tolerances`() {
        // On [0, 1]: x_1 = n^{-6} > 4e-15 iff n <= 250 (250^{-6} = 4.10e-15, 251^{-6} = 3.999e-15).
        Grid.power(250, r = 6.0)
        assertFailsWith<IllegalArgumentException> { Grid.power(251, r = 6.0) }
        assertFailsWith<IllegalArgumentException> { Grid.power(256, r = 6.0) }
        val e = assertFailsWith<IllegalArgumentException> { Grid.power(512, r = 6.0) }
        assertTrue(e.message!!.contains("breakpointInclusionEps"), "message: ${e.message}")
        // The tolerance scales with b - a, so a longer interval does not move the threshold.
        assertFailsWith<IllegalArgumentException> { Grid.power(256, 0.0, 1e3, 6.0) }
        Grid.symmetricPower(64, r = 6.0)
        assertFailsWith<IllegalArgumentException> { Grid.symmetricPower(512, r = 6.0) }
        val powerFirstReject = (1..600).first { n -> runCatching { Grid.power(n, r = 6.0) }.isFailure }
        val symFirstReject = (2..600 step 2).first { n -> runCatching { Grid.symmetricPower(n, r = 6.0) }.isFailure }
        record("guard r=6 on [0,1]: power throws from n=$powerFirstReject; symmetricPower throws from n=$symFirstReject")
        assertEquals(251, powerFirstReject)
    }

    @Test
    fun `invalid arguments are rejected`() {
        assertFailsWith<IllegalArgumentException> { Grid.power(0, r = 2.0) }
        assertFailsWith<IllegalArgumentException> { Grid.power(-1, r = 2.0) }
        assertFailsWith<IllegalArgumentException> { Grid.power(8, r = 0.5) }
        assertFailsWith<IllegalArgumentException> { Grid.power(8, r = Double.NaN) }
        assertFailsWith<IllegalArgumentException> { Grid.power(8, r = Double.POSITIVE_INFINITY) }
        assertFailsWith<IllegalArgumentException> { Grid.power(8, 1.0, 0.0, 2.0) }
        assertFailsWith<IllegalArgumentException> { Grid.power(8, Double.NaN, 1.0, 2.0) }
        assertFailsWith<IllegalArgumentException> { Grid.symmetricPower(0, r = 2.0) }
        assertFailsWith<IllegalArgumentException> { Grid.symmetricPower(7, r = 2.0) }
        assertFailsWith<IllegalArgumentException> { Grid.symmetricPower(9, r = 2.0) }
        assertFailsWith<IllegalArgumentException> { Grid.symmetricPower(8, r = 0.99) }
        assertFailsWith<IllegalArgumentException> { Grid.symmetricPower(8, r = Double.NaN) }
        assertFailsWith<IllegalArgumentException> { Grid.symmetricPower(8, 0.0, Double.POSITIVE_INFINITY, 2.0) }
    }

    @Test
    fun `minimal spline bases on power-graded grids are biorthogonal to their projection functionals`() {
        val cases = listOf(
            "B, power(64, r=4)" to MinimalSplineBasis(GeneratingSystem.B, Grid.power(64, r = 4.0)),
            "B, symmetricPower(64, r=3)" to MinimalSplineBasis(GeneratingSystem.B, Grid.symmetricPower(64, r = 3.0)),
            "pow(0.5), power(64, r=4)" to MinimalSplineBasis(
                GeneratingSystem.reparametrized(Reparametrization.power(0.5)), Grid.power(64, r = 4.0),
            ),
        )
        for ((name, basis) in cases) {
            val n = basis.grid.n
            val theta = ProjFunctionals(basis)
            var bio = 0.0
            for (i in -2..n - 1) for (j in -2..n - 1) {
                val v = theta.chi(i).apply { t -> basis.omega(j, t) }
                bio = max(bio, abs(v - if (i == j) 1.0 else 0.0))
            }
            record("biorthogonality $name: $bio")
            assertTrue(bio <= 1e-12, "biorthogonality $name: $bio")
        }
    }
}
