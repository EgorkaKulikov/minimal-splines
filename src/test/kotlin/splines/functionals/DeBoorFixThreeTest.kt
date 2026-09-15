package splines.functionals

import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertTrue

/**
 * Tests of the three de Boor--Fix functionals xi^<0>, xi^<1>, xi^<2>
 * The sources of the formulas are listed in docs/REFERENCES.md, section 2.
 */
@Tag("fast")
class DeBoorFixThreeTest {
    private val grids = listOf(Grid.uniform(8), Grid.quasiUniform(8))
    private val systems = listOf(GeneratingSystem.B, GeneratingSystem.H, GeneratingSystem.T)

    /** (a) Biorthogonality xi^<r>_i(omega_j)=delta_ij for r=0,1,2, all of B/H/T, ALL i,j (boundary ones included). */
    @Test fun biorthogonalityAllRAllBasesIncludingBoundary() {
        for (grid in grids) for (sys in systems) {
            val basis = MinimalSplineBasis(sys, grid)
            val n = basis.n
            for (r in 0..2) {
                val xi = DeBoorFixFunctionals(basis, r)
                for (i in -2 until n) for (j in -2 until n) {
                    val v = xi.chi(i).apply(
                        { t -> basis.omega(j, t) },
                        { t -> basis.omegaDeriv(j, t) },
                        { t -> basis.omegaDeriv2(j, t) },
                    )
                    val expected = if (i == j) 1.0 else 0.0
                    assertTrue(
                        kotlin.math.abs(v - expected) < 1e-8,
                        "xi<$r>_$i(omega_$j)=$v (sys=${sys.name}, expected=$expected)",
                    )
                }
            }
        }
    }

    /** (b) Exactness of the projector P_xi on span{1,rho,sigma} for every r and basis. */
    @Test fun projectorExactOnGeneratingSpan() {
        val ts = (0..100).map { it / 100.0 }
        for (grid in grids) for (sys in systems) {
            val basis = MinimalSplineBasis(sys, grid)
            val members = listOf(
                Triple<(Double) -> Double, (Double) -> Double, (Double) -> Double>(
                    { 1.0 }, { 0.0 }, { 0.0 },
                ) to "1",
                Triple(sys.rho, sys.rhoD, sys.rhoDD) to "rho",
                Triple(sys.sigma, sys.sigmaD, sys.sigmaDD) to "sigma",
            )
            for (r in 0..2) {
                val xi = DeBoorFixFunctionals(basis, r)
                for ((triple, label) in members) {
                    val (g, gD, gDD) = triple
                    val pc = xi.projectorCoeffs(g, gD, gDD)
                    for (ti in ts) {
                        val t = grid.a + (grid.b - grid.a) * ti
                        val err = kotlin.math.abs(g(t) - basis.evalSpline(pc, t))
                        assertTrue(
                            err < 1e-8,
                            "P_xi<$r> $label at t=$t err=$err (sys=${sys.name})",
                        )
                    }
                }
            }
        }
    }
}
