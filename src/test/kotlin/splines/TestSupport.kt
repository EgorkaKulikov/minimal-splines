package splines

import splines.metrics.DEFAULT_CONTROL_REFINEMENT
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.exp
import kotlin.math.max
import kotlin.math.sin

/**
 * Maximum of |exact − eval| over a control grid of `refinement·n + 1` points, excluding a
 * boundary layer of width `layer·h` on each side (h is the maximal grid step).
 *
 * Needed for the ξ̃ family (`DiscreteDeBoorFixFunctionals`): at multiple boundary nodes the central
 * difference degenerates into a one-sided one, and the quasi-projector does not reproduce span φ in
 * the boundary layer (defect ≈ h²). This is a property of the construction, not an implementation
 * defect, hence exactness on span φ and order 3 for ξ̃ are checked in the interior region.
 */
internal fun interiorErrorEh(
    exact: (Double) -> Double,
    eval: (Double) -> Double,
    grid: Grid,
    layer: Int = 3,
    refinement: Int = DEFAULT_CONTROL_REFINEMENT,
): Double {
    val m = refinement * grid.n
    val lo = grid.a + layer * grid.h
    val hi = grid.b - layer * grid.h
    var mx = 0.0
    for (i in 0..m) {
        val t = grid.a + (grid.b - grid.a) * i / m
        if (t < lo || t > hi) continue
        mx = max(mx, abs(exact(t) - eval(t)))
    }
    return mx
}

/** Test function f, f', f'' as a triple of lambdas. */
internal data class TestFunction(
    val f: (Double) -> Double,
    val fD: (Double) -> Double,
    val fDD: (Double) -> Double,
)

/**
 * f(t) = exp(sin 3t), which lies in span φ of none of the systems B, H, T;
 * f' = 3 cos 3t · f, f'' = (−9 sin 3t + 9 cos² 3t) · f.
 */
internal fun testFunction(): TestFunction = TestFunction(
    f = { t -> exp(sin(3.0 * t)) },
    fD = { t -> 3.0 * cos(3.0 * t) * exp(sin(3.0 * t)) },
    fDD = { t ->
        val s = sin(3.0 * t)
        val c = cos(3.0 * t)
        (-9.0 * s + 9.0 * c * c) * exp(s)
    },
)
