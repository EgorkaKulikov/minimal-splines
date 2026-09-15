package splines

/** Explicit formulas for quadratic B- and H-splines, for cross-checking the general construction. */
object ReferenceSplines {
    /** Classical quadratic B-spline omega^B_j(t). */
    fun omegaB(grid: Grid, j: Int, t: Double): Double {
        val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2); val xj3 = grid.x(j + 3)
        if (t < xj || t > xj3) return 0.0
        return when {
            t < xj1 -> (t - xj) * (t - xj) / ((xj1 - xj) * (xj2 - xj))
            t < xj2 -> (1.0 / (xj1 - xj)) * (
                (t - xj) * (t - xj) / (xj2 - xj)
                    - (t - xj1) * (t - xj1) * (xj3 - xj) / ((xj2 - xj1) * (xj3 - xj1))
                )
            else -> (t - xj3) * (t - xj3) / ((xj3 - xj1) * (xj3 - xj2))
        }
    }

    /** Derivative of the golden reference B-spline. */
    fun omegaBDeriv(grid: Grid, j: Int, t: Double): Double {
        val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2); val xj3 = grid.x(j + 3)
        if (t < xj || t > xj3) return 0.0
        return when {
            t < xj1 -> 2.0 * (t - xj) / ((xj1 - xj) * (xj2 - xj))
            t < xj2 -> (1.0 / (xj1 - xj)) * (
                2.0 * (t - xj) / (xj2 - xj)
                    - 2.0 * (t - xj1) * (xj3 - xj) / ((xj2 - xj1) * (xj3 - xj1))
                )
            else -> 2.0 * (t - xj3) / ((xj3 - xj1) * (xj3 - xj2))
        }
    }

    /** Hyperbolic minimal spline omega^H_j(t). */
    fun omegaH(grid: Grid, j: Int, t: Double): Double {
        val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2); val xj3 = grid.x(j + 3)
        if (t < xj || t > xj3) return 0.0
        fun sh(x: Double) = Math.sinh(x)
        fun ch(x: Double) = Math.cosh(x)
        val cosCenter = ch((xj2 - xj1) / 2.0)
        return when {
            t < xj1 -> cosCenter * sh((t - xj) / 2.0) * sh((t - xj) / 2.0) /
                (sh((xj1 - xj) / 2.0) * sh((xj2 - xj) / 2.0))
            t < xj2 -> (cosCenter / sh((xj1 - xj) / 2.0)) * (
                sh((t - xj) / 2.0) * sh((t - xj) / 2.0) / sh((xj2 - xj) / 2.0)
                    - sh((xj3 - xj) / 2.0) * sh((t - xj1) / 2.0) * sh((t - xj1) / 2.0) /
                    (sh((xj3 - xj1) / 2.0) * sh((xj2 - xj1) / 2.0))
                )
            else -> cosCenter * sh((xj3 - t) / 2.0) * sh((xj3 - t) / 2.0) /
                (sh((xj3 - xj1) / 2.0) * sh((xj3 - xj2) / 2.0))
        }
    }
}
