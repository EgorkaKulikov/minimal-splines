package splines.metrics

import kotlin.math.abs

import splines.Grid

// ============================================================================
// Error metrics
// ============================================================================

/**
 * Default refinement factor of the control grid in [errorEh]: `100 n` control intervals,
 * that is, `100 n + 1` control points.
 */
public const val DEFAULT_CONTROL_REFINEMENT: Int = 100

/**
 * Uniform error norm `E_h = max |exact(t) - eval(t)|` over `refinement * n + 1` points of the
 * uniform control grid on the interval `[a, b]` of the grid [grid].
 *
 * The maximum over a finite set of points does not exceed the exact maximum; the default value of
 * [refinement] is checked to be sufficient by recomputing on a finer control grid.
 *
 * @param exact exact function.
 * @param eval approximation (a spline, for example).
 * @param grid grid defining the interval `[a, b]` and the number of intervals `n`.
 * @param refinement number of control intervals per grid interval.
 * @throws IllegalArgumentException if `grid.n < 1` or `refinement < 1`: the number of control
 *   intervals `m = refinement * n` is used as a divisor when computing the control points.
 */
public fun errorEh(
    exact: (Double) -> Double,
    eval: (Double) -> Double,
    grid: Grid,
    refinement: Int = DEFAULT_CONTROL_REFINEMENT,
): Double {
    require(grid.n >= 1) { "errorEh: grid.n must be at least 1 (otherwise m = refinement*n = 0), got n=${grid.n}" }
    require(refinement >= 1) { "errorEh: refinement must be at least 1, got $refinement" }
    val m = refinement * grid.n
    var e = 0.0
    for (i in 0..m) {
        // The product (b - a)·i/m at i = m may round to the right of b by one unit in the last
        // place (for example, b = 0.013, m = 130); the spline is not defined outside [a, b],
        // so the control point is clamped to the grid interval. For points inside the
        // interval the clamping does not change a single bit.
        val t = (grid.a + (grid.b - grid.a) * i / m).coerceIn(grid.a, grid.b)
        e = maxOf(e, abs(exact(t) - eval(t)))
    }
    return e
}
