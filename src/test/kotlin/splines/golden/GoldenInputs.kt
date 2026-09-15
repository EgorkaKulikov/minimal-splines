package splines.golden

import splines.GeneratingSystem
import splines.Grid
import kotlin.math.PI
import kotlin.math.cos
import kotlin.math.exp
import kotlin.math.sin

/**
 * Deterministic inputs of the golden tests. No randomness: all grids, points,
 * functions and coefficients are given by explicit formulas, so that the golden reference and the
 * check are built from the same values on any platform.
 */
object GoldenInputs {
    val NS: List<Int> = listOf(1, 2, 4, 8, 16, 32)

    /** All grids in a deterministic order; the key has the form `type-n` (plus two special ones). */
    val grids: LinkedHashMap<String, Grid> = linkedMapOf<String, Grid>().apply {
        for (n in NS) {
            put("uniform-$n", Grid.uniform(n, 0.0, 1.0))
            put("quasiUniform-$n", Grid.quasiUniform(n, 0.0, 1.0))
            if (n >= 2) {
                put("geometric-$n", Grid.geometric(n, 0.0, 1.0))
                put("graded-$n", Grid.graded(n, 0.0, 1.0))
            }
        }
        // Interval of length 4 > pi: the T system cannot be built here (the Wronskian changes sign).
        put("uniform-8-wide", Grid.uniform(8, -1.5, 2.5))
        // Interval of length 2 < pi: an extra grid for the T system.
        put("uniform-8-T", Grid.uniform(8, 0.0, 2.0))
    }

    val systems: LinkedHashMap<String, GeneratingSystem> = linkedMapOf(
        "B" to GeneratingSystem.B,
        "H" to GeneratingSystem.H,
        "T" to GeneratingSystem.T,
    )

    /** T is admissible only on intervals of length strictly less than pi. */
    fun allowed(sysName: String, grid: Grid): Boolean = sysName != "T" || (grid.b - grid.a) < PI

    /** Pairs (system, grid) for which the basis is defined. */
    fun basisCases(maxN: Int = Int.MAX_VALUE): List<Triple<String, String, Grid>> = buildList {
        for ((sName, _) in systems) for ((gName, g) in grids) {
            if (g.n <= maxN && allowed(sName, g)) add(Triple(sName, gName, g))
        }
    }

    /** 7 points strictly inside every grid interval: x_k + (x_{k+1} - x_k)(i + 0.5)/7. */
    fun controlPoints(grid: Grid): DoubleArray {
        val out = DoubleArray(7 * grid.n)
        var p = 0
        for (k in 0 until grid.n) {
            val xk = grid.x(k); val xk1 = grid.x(k + 1)
            for (i in 0 until 7) out[p++] = xk + (xk1 - xk) * (i + 0.5) / 7.0
        }
        return out
    }

    class TestFunction(
        val name: String,
        val f: (Double) -> Double,
        val fD: (Double) -> Double,
        val fDD: (Double) -> Double,
    )

    val functions: List<TestFunction> = listOf(
        // f1 = exp(sin 3t); f1' = 3 cos(3t) f1; f1'' = (9 cos^2(3t) - 9 sin(3t)) f1
        TestFunction(
            "f1",
            { t -> exp(sin(3.0 * t)) },
            { t -> 3.0 * cos(3.0 * t) * exp(sin(3.0 * t)) },
            { t -> (9.0 * cos(3.0 * t) * cos(3.0 * t) - 9.0 * sin(3.0 * t)) * exp(sin(3.0 * t)) },
        ),
        // f2 = 1/(1+t^2); f2' = -2t/(1+t^2)^2; f2'' = (6t^2 - 2)/(1+t^2)^3
        TestFunction(
            "f2",
            { t -> 1.0 / (1.0 + t * t) },
            { t -> val q = 1.0 + t * t; -2.0 * t / (q * q) },
            { t -> val q = 1.0 + t * t; (6.0 * t * t - 2.0) / (q * q * q) },
        ),
        // f3 = t^3 - 2t + 0.5; f3' = 3t^2 - 2; f3'' = 6t
        TestFunction(
            "f3",
            { t -> t * t * t - 2.0 * t + 0.5 },
            { t -> 3.0 * t * t - 2.0 },
            { t -> 6.0 * t },
        ),
        // f4 = cos(2t) e^{-t}; f4' = -e^{-t}(2 sin 2t + cos 2t); f4'' = e^{-t}(4 sin 2t - 3 cos 2t)
        TestFunction(
            "f4",
            { t -> cos(2.0 * t) * exp(-t) },
            { t -> -exp(-t) * (2.0 * sin(2.0 * t) + cos(2.0 * t)) },
            { t -> exp(-t) * (4.0 * sin(2.0 * t) - 3.0 * cos(2.0 * t)) },
        ),
    )

    /** Spline coefficients of length n+2: c_j = sin(1.7 j + 0.3). */
    fun coeffs(n: Int): DoubleArray = DoubleArray(n + 2) { sin(1.7 * it + 0.3) }

    /** 20 points on [-0.5, 2] for phi/phiD/phiDD/wronskian. */
    val phiPoints: DoubleArray = DoubleArray(20) { -0.5 + 2.5 * (it + 0.5) / 20.0 }
}
