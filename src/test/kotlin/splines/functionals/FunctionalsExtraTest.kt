package splines.functionals

import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import org.junit.jupiter.api.Tag
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Tests of the quasi-projection families (mu, lambda) and of the shared infrastructure
 * FunctionalFamily/ValueFunctional/DerivFunctional. The main invariant of
 * quasi-interpolants is exactness on span{1,rho,sigma}: P_chi g = g.
 */
@Tag("fast")
class FunctionalsExtraTest {
    private val grid = Grid.uniform(8)
    private val basis = MinimalSplineBasis(GeneratingSystem.B, grid)

    /** Reproduction of g by the spline P_chi g = sum chi_j(g) omega_j at a set of points. */
    private fun assertReproduces(fam: FunctionalFamily, g: (Double) -> Double, tol: Double) {
        val c = fam.projectorCoeffs(g)
        for (t in listOf(0.07, 0.23, 0.5, 0.71, 0.93)) {
            assertEquals(g(t), basis.evalSpline(c, t), tol, "${fam.name}: P g != g at t=$t")
        }
    }

    /** mu (averaging): exact on 1, rho=t, sigma=t^2 for phi^B (quasi-interpolant of the span). */
    @Test fun averagingExactOnSpanB() {
        val mu = AveragingFunctionals(basis)
        assertTrue(!mu.isProjector && !mu.usesDerivative)
        assertReproduces(mu, { 1.0 }, 1e-9)
        assertReproduces(mu, { t -> t }, 1e-8)
        assertReproduces(mu, { t -> t * t }, 1e-8)
    }

    /** mu on the hyperbolic system H: exactness on span{1,sinh,cosh}. */
    @Test fun averagingExactOnSpanH() {
        val basisH = MinimalSplineBasis(GeneratingSystem.H, grid)
        val mu = AveragingFunctionals(basisH)
        val c = mu.projectorCoeffs({ t -> Math.sinh(t) })
        for (t in listOf(0.2, 0.6, 0.85)) {
            assertEquals(Math.sinh(t), basisH.evalSpline(c, t), 1e-7, "mu H sinh at $t")
        }
    }

    /** lambda (three-point): exact on span{1,rho,sigma} for phi^B. */
    @Test fun threePointExactOnSpanB() {
        val lam = ThreePointFunctionals(basis)
        assertTrue(!lam.isProjector && !lam.usesDerivative)
        assertReproduces(lam, { 1.0 }, 1e-9)
        assertReproduces(lam, { t -> t }, 1e-8)
        assertReproduces(lam, { t -> t * t }, 1e-8)
    }

    /** The boundary mu_{-2}, mu_{n-1} are pure values u(x_0), u(x_n). */
    @Test fun averagingBoundaryAreEndpointValues() {
        val mu = AveragingFunctionals(basis)
        assertEquals(grid.x(0) * 3.0, mu.chi(-2).apply { t -> 3.0 * t }, 1e-12)
        assertEquals(grid.x(grid.n) * 3.0, mu.chi(grid.n - 1).apply { t -> 3.0 * t }, 1e-12)
        assertEquals(1.0, mu.chi(-2).absSum(), 1e-12)
    }

    /** The boundary lambda_{-2}, lambda_{n-1} are pure values at the endpoints. */
    @Test fun threePointBoundaryAreEndpointValues() {
        val lam = ThreePointFunctionals(basis)
        assertEquals(5.0, lam.chi(-2).apply { _ -> 5.0 }, 1e-12)
        assertEquals(5.0, lam.chi(grid.n - 1).apply { _ -> 5.0 }, 1e-12)
    }

    /** cChi() = max_j sum|coeff| is a positive finite stability constant. */
    @Test fun cChiPositiveFinite() {
        for (fam in listOf(AveragingFunctionals(basis), ThreePointFunctionals(basis), ProjFunctionals(basis))) {
            val c = fam.cChi()
            assertTrue(c >= 1.0 && c.isFinite(), "${fam.name}: cChi=$c")
        }
    }

    /** projectorCoeffs returns a vector of length n+2. */
    @Test fun projectorCoeffsLength() {
        val mu = AveragingFunctionals(basis)
        assertEquals(grid.n + 2, mu.projectorCoeffs({ t -> t }).size)
    }

    /** DerivFunctional: apply = f(node)+cD*f'(node); absSum = 1+|cD|. */
    @Test fun derivFunctionalApplyAndAbsSum() {
        val df = DerivFunctional(0.5, -2.0)
        // f(t)=t^2 -> f(0.5)=0.25, f'(0.5)=1.0 -> 0.25 + (-2)*1 = -1.75
        assertEquals(-1.75, df.apply({ t -> t * t }, { t -> 2 * t }), 1e-12)
        assertEquals(3.0, df.absSum(), 1e-12) // 1+|-2|
    }

    /** ValueFunctional.absSum = sum of the absolute values of the coefficients. */
    @Test fun valueFunctionalAbsSum() {
        val vf = ValueFunctional(doubleArrayOf(0.0, 1.0), doubleArrayOf(-3.0, 4.0))
        assertEquals(7.0, vf.absSum(), 1e-12)
    }

    /** xi family (de Boor–Fix): apply without a derivative ignores the derivative term at the boundaries. */
    @Test fun deBoorFixUsesDerivativeFlag() {
        val xi = DeBoorFixFunctionals(basis)
        assertTrue(xi.isProjector && xi.usesDerivative)
        // the boundary xi_{-2} = u(x_0): a pure value, no derivative needed
        assertEquals(2.0 * grid.x(0), xi.chi(-2).apply { t -> 2.0 * t }, 1e-12)
    }

    /** closedFormInternal agrees with the biorthogonal theta on omega_i (health check). */
    @Test fun closedFormMatchesBuiltTheta() {
        val theta = ProjFunctionals(basis)
        val j = 2
        val closed = theta.closedFormInternal(j)
        // theta_j(omega_i)=delta_ij holds through the closed form as well
        for (i in (j - 2)..(j + 2)) {
            val v = closed.apply { t -> basis.omega(i, t) }
            assertEquals(if (i == j) 1.0 else 0.0, v, 1e-6, "closedForm theta_$j(omega_$i)")
        }
    }
}
