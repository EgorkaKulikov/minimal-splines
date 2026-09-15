package splines.property

import net.jqwik.api.Arbitraries
import net.jqwik.api.Arbitrary
import net.jqwik.api.Combinators
import net.jqwik.api.RandomDistribution
import splines.GeneratingSystem
import splines.Grid

/** A random grid together with its system: for T the interval length is bounded (the Wronskian keeps its sign). */
data class GridCase(val sys: GeneratingSystem, val grid: Grid) {
    override fun toString(): String =
        "GridCase(sys=${sys.name}, n=${grid.n}, a=${grid.a}, b=${grid.b}, h=${grid.h})"
}

object Generators {
    private const val MAX_LENGTH_TRIG = 3.0

    fun systems(): Arbitrary<GeneratingSystem> =
        Arbitraries.of(GeneratingSystem.B, GeneratingSystem.H, GeneratingSystem.T)

    /**
     * Interior nodes in [0, 1] built from n positive steps whose neighbouring ratio is at most
     * 1000: the relative gap is automatically >= 1/(1000·n) > 1e-3/n, with no rejection (a gap filter
     * on jqwik's biased generator produced TooManyFilterMisses).
     */
    private fun unitInterior(n: Int): Arbitrary<DoubleArray> =
        Arbitraries.doubles().between(1.0, 1000.0).ofScale(3)
            .withDistribution(RandomDistribution.uniform())
            .array(DoubleArray::class.java).ofSize(n)
            .map { gaps ->
                val cum = DoubleArray(n + 1)
                for (i in 0 until n) cum[i + 1] = cum[i] + gaps[i]
                DoubleArray(n - 1) { cum[it + 1] / cum[n] }
            }

    fun grids(sys: GeneratingSystem): Arbitrary<Grid> {
        val maxLen = if (sys === GeneratingSystem.T) MAX_LENGTH_TRIG else 10.0
        return Arbitraries.integers().between(3, 40).flatMap { n ->
            Combinators.combine(
                Arbitraries.doubles().between(-10.0, 10.0),
                Arbitraries.doubles().between(1e-3, maxLen).ofScale(3),
                unitInterior(n),
            ).`as` { a, len, unit ->
                val pts = DoubleArray(n + 1)
                pts[0] = a
                pts[n] = a + len
                for (i in 1 until n) pts[i] = a + len * unit[i - 1]
                Grid(n, pts)
            }
        }
    }

    fun gridCases(): Arbitrary<GridCase> = systems().flatMap { sys -> grids(sys).map { GridCase(sys, it) } }

    fun coefficients(): Arbitrary<DoubleArray> =
        Arbitraries.doubles().between(-5.0, 5.0).array(DoubleArray::class.java).ofSize(3)
}
