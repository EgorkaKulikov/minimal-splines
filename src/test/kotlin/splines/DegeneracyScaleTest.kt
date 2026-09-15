package splines

import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Scale invariance of the degeneracy criteria: the basis is built on intervals of length 1e-6
 * and 1e6. The criteria (the relative threshold [DEGENERACY_RELATIVE_EPS] for scalar denominators
 * and the condition number of the approximation relation matrix in local coordinates) do not
 * depend on the scale of the interval: a valid grid on a small interval is not rejected, and on a
 * large interval rounding noise is not mistaken for a significant quantity.
 */
@Tag("fast")
class DegeneracyScaleTest {

    /**
     * Small scale: the minimal spline basis on the interval `[0, 1e-6]` for the polynomial system B
     * is built, and the partition of unity `sum_j omega_j(t) = 1` holds (the first component of phi
     * equals 1). An absolute degeneracy threshold would reject such a grid, since
     * `det(M_k) ~ h^3 ~ 2e-21` in global coordinates.
     */
    @Test
    fun splineBasisBuildsOnTinyInterval() {
        val grid = Grid.uniform(8, 0.0, 1e-6)
        val basis = MinimalSplineBasis(GeneratingSystem.B, grid)
        for (k in 0 until grid.n) {
            val t = 0.5 * (grid.x(k) + grid.x(k + 1))
            var sum = 0.0
            for (j in -2..grid.n - 1) sum += basis.omega(j, t)
            assertEquals(1.0, sum, 1e-9, "partition of unity at t=$t")
        }
    }

    /**
     * Non-polynomial systems (H, T) on an interval of length `1e-6`: either the construction succeeds,
     * and then the partition of unity holds, or it is rejected with a diagnostic that names the reason
     * (singularity or ill conditioning). Returning unreliable values is excluded.
     */
    @Test
    fun nonPolynomialSystemOnTinyInterval() {
        for (sys in listOf(GeneratingSystem.H, GeneratingSystem.T)) {
            val grid = Grid.uniform(8, 0.0, 1e-6)
            val basis = try { MinimalSplineBasis(sys, grid) } catch (ex: IllegalArgumentException) {
                println("DegeneracyScaleTest: ${sys.name} on [0, 1e-6] rejected: ${ex.message}")
                assertTrue(
                    ex.message!!.contains("ill-conditioned") || ex.message!!.contains("numerically singular"),
                    "the diagnostic must name the reason for the rejection: ${ex.message}"
                )
                continue
            }
            println("DegeneracyScaleTest: ${sys.name} on [0, 1e-6] built")
            val ones = DoubleArray(grid.n + 2) { 1.0 }
            for (i in 1..20) {
                val t = 1e-6 * i / 21.0
                assertEquals(1.0, basis.evalSpline(ones, t), 1e-8, "partition of unity for ${sys.name} at t=$t")
            }
        }
    }

    /**
     * Large scale: on `[0, 1e6]` the construction succeeds, since the grid is non-degenerate.
     *
     * The polynomial and trigonometric systems are used: for the hyperbolic one `cosh(1e6)`
     * overflows double, which is a limitation of the generating system and not of the degeneracy
     * criterion.
     */
    @Test
    fun splineBasisBuildsOnHugeInterval() {
        for (sys in listOf(GeneratingSystem.B, GeneratingSystem.T)) {
            val grid = Grid.uniform(8, 0.0, 1e6)
            val basis = MinimalSplineBasis(sys, grid)
            for (k in 0 until grid.n) {
                val t = 0.5 * (grid.x(k) + grid.x(k + 1))
                var sum = 0.0
                for (j in -2..grid.n - 1) sum += basis.omega(j, t)
                assertEquals(1.0, sum, 1e-6, "partition of unity for ${sys.name} at t=$t")
            }
        }
    }
}
