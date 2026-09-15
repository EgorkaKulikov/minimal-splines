package splines.functionals

import kotlin.math.abs

import numerics.DenseMatrix
import numerics.LinearAlgebra
import numerics.NumericsContext
import splines.DEGENERACY_RELATIVE_EPS
import splines.Grid
import splines.MinimalSplineBasis
import splines.cancellationScale
import splines.isSignificant

// ============================================================================
// Families of (quasi-)projection functionals: theta (projection), xi (de Boor–Fix,
// value and derivatives), xitilde (discretized de Boor–Fix), mu (averaging),
// lambda (three-point). The sources of the formulas are given in the section headers.
// ============================================================================

/**
 * Approximation functional chi_j. The theta, mu and lambda families use only the values of f;
 * the xi family also uses the derivatives f' and f'', so the interface accepts the function
 * together with its derivatives, while the families without derivatives ignore them.
 */
public interface ApproxFunctional {
    /**
     * Value chi_j(f) for the function f, its first derivative fD and second derivative fDD.
     * The theta, mu and lambda families use only f; xi^<1>, xi^<2> also use fD; xi^<0> also uses fDD.
     * By default fDD is zero, which is sufficient for the families that do not use the second derivative.
     */
    public fun apply(f: (Double) -> Double, fD: (Double) -> Double, fDD: (Double) -> Double = { 0.0 }): Double

    /**
     * Sum of the absolute values of the functional coefficients in its representation, without
     * normalization by the grid step.
     *
     * For value functionals (`chi_j(f) = sum_k c_k f(t_k)`; the theta, mu, lambda and xitilde
     * families) this quantity coincides with the norm of the functional with respect to
     * perturbations of the input values: a data perturbation of size eps in the maximum norm
     * produces an error of at most `absSum() * eps`.
     * For functionals involving derivatives ([DerivFunctional], [SecondDerivFunctional]) the
     * quantity is not such a norm (see the KDoc of those classes) and serves only as a diagnostic
     * of the representation.
     */
    public fun absSum(): Double
}

/** Convenience wrapper: chi_j(f) without explicit derivatives (derivatives = 0). */
public fun ApproxFunctional.apply(f: (Double) -> Double): Double = apply(f, { 0.0 }, { 0.0 })

/**
 * Family of approximation functionals {chi_j}_{j=-2}^{n-1} and the (quasi-)projector
 * P_chi g = sum_j chi_j(g) omega_j. Common interface of the theta, xi, xitilde, mu and lambda families.
 *
 * @property basis minimal spline basis the functionals belong to.
 * @property name family name (theta, xi, xi<r>, xitilde, xitilde<r>, mu, lambda).
 * @property isProjector `true` for projectors (theta, xi): biorthogonality
 *   chi_i(omega_j) = delta_ij holds and P_chi^2 = P_chi; `false` for quasi-interpolants (xitilde, mu, lambda).
 * @property usesDerivative `true` for the xi family, which uses the derivative of its argument.
 */
public abstract class FunctionalFamily(
    public val basis: MinimalSplineBasis,
    public val name: String,
    /**
     * Numerical computation context: the BLAS/LAPACK implementation the theta, mu and lambda
     * families use to solve the 3×3 and 5×5 linear systems in the constructor. All families accept
     * the parameter, including those that do not use linear algebra (xi, xitilde), which allows
     * the context of any family to be checked against the context of the calling code uniformly.
     */
    public val ctx: NumericsContext = NumericsContext.default(),
) {
    /** Grid of the basis. */
    public val grid: Grid = basis.grid

    /** Number of grid intervals; the family consists of n + 2 functionals. */
    public val n: Int = grid.n
    public abstract val isProjector: Boolean
    public abstract val usesDerivative: Boolean

    /** `true` for families that use the second derivative of their argument (xi^<0>). */
    public open val usesSecondDerivative: Boolean = false

    /** Functional chi_j, j = -2..n-1. */
    public abstract fun chi(j: Int): ApproxFunctional

    /** Coefficients of the projection P_chi g = sum chi_j(g) omega_j: the vector (chi_j(g)) of size n+2. */
    public fun projectorCoeffs(
        g: (Double) -> Double,
        gD: (Double) -> Double = { 0.0 },
        gDD: (Double) -> Double = { 0.0 },
    ): DoubleArray = DoubleArray(n + 2) { chi(it - 2).apply(g, gD, gDD) }

    /**
     * Maximum of [ApproxFunctional.absSum] over all j = -2..n-1.
     *
     * For families of value functionals (`usesDerivative == false`) this is the constant `C_chi` —
     * a bound on the amplification of an input data perturbation by the (quasi-)projector P_chi in
     * the maximum norm. For the xi family (`usesDerivative == true`) the quantity is not such a
     * bound: the coefficients at the derivatives decay with the grid step, whereas the
     * amplification of a perturbation by numerical differentiation grows (see the KDoc of
     * [DerivFunctional]).
     */
    public fun cChi(): Double = (-2..n - 1).maxOf { chi(it).absSum() }
}

// ----------------------------------------------------------------------------
// theta — projection functionals
// Source: Kulikov, Makarov (Zapiski Nauchnykh Seminarov POMI, 2025, vol. 542,
// p. 126–143). See docs/REFERENCES.md, section 2.
// ----------------------------------------------------------------------------

/**
 * Value functional: a linear combination of the values of f at the points [nodes] with the
 * coefficients [coeffs].
 *
 * The arrays are not copied; modifying their contents breaks the invariants of the functional.
 *
 * @property nodes support points.
 * @property coeffs coefficients at the values at [nodes].
 * @throws IllegalArgumentException if the array lengths differ.
 */
public class ValueFunctional(public val nodes: DoubleArray, public val coeffs: DoubleArray) : ApproxFunctional {
    init { require(nodes.size == coeffs.size) }
    override fun apply(f: (Double) -> Double, fD: (Double) -> Double, fDD: (Double) -> Double): Double {
        var s = 0.0
        for (k in nodes.indices) s += coeffs[k] * f(nodes[k])
        return s
    }
    /**
     * Sum of the absolute values of [coeffs]. For a value functional this is exactly the norm of
     * the functional as the amplification factor of a perturbation of the input values:
     * `|chi(f + e) - chi(f)| <= absSum() * max|e|`.
     */
    override fun absSum(): Double = coeffs.fold(0.0) { acc, v -> acc + abs(v) }
}

/**
 * Family of projection functionals theta_j.
 *
 * The internal and boundary functionals are built by local biorthogonalization: the system
 * theta_j(omega_i) = delta_ij is solved at the grid nodes and the interval midpoints. This is a
 * stable equivalent representation of the closed-form formula from the source; the agreement of
 * the two representations is checked by a test (see [closedFormInternal]).
 *
 * The family is a projector: biorthogonality theta_i(omega_j) = delta_ij holds, hence
 * P_theta^2 = P_theta. The boundary functionals (j = -2 and j = n-1) are the values f(x_0) and
 * f(x_n) according to the definition in the source.
 *
 * @param ctx numerical computation context for the linear systems of the local biorthogonalization.
 */
public class ProjFunctionals(
    basis: MinimalSplineBasis,
    ctx: NumericsContext = NumericsContext.default(),
) : FunctionalFamily(basis, "theta", ctx) {
    override val isProjector: Boolean = true
    override val usesDerivative: Boolean = false
    private val funcs: Array<ApproxFunctional> = Array(n + 2) { buildTheta(it - 2) }
    override fun chi(j: Int): ApproxFunctional = funcs[j + 2]

    private fun mid(p: Double, q: Double) = 0.5 * (p + q)

    private fun buildTheta(j: Int): ApproxFunctional {
        val x0 = grid.x(0); val x1 = grid.x(1)
        val xnm1 = grid.x(n - 1); val xn = grid.x(n)
        return when (j) {
            -2 -> ValueFunctional(doubleArrayOf(x0), doubleArrayOf(1.0))
            -1 -> localFunctional(j = -1, points = doubleArrayOf(x0, mid(x0, x1), x1), indices = intArrayOf(-2, -1, 0))
            n - 1 -> ValueFunctional(doubleArrayOf(xn), doubleArrayOf(1.0))
            n - 2 -> localFunctional(j = n - 2, points = doubleArrayOf(xnm1, mid(xnm1, xn), xn), indices = intArrayOf(n - 3, n - 2, n - 1))
            else -> {
                val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2); val xj3 = grid.x(j + 3)
                localFunctional(
                    j = j,
                    points = doubleArrayOf(xj, mid(xj, xj1), mid(xj1, xj2), mid(xj2, xj3), xj3),
                    indices = intArrayOf(j - 2, j - 1, j, j + 1, j + 2),
                )
            }
        }
    }

    /** Local biorthogonalization: coeff such that sum_p coeff_p omega_i(points_p)=delta_ij. */
    private fun localFunctional(j: Int, points: DoubleArray, indices: IntArray): ValueFunctional {
        val m = points.size
        val matrix = DenseMatrix.build(m, m) { r, c -> basis.omega(indices[r], points[c]) }
        val rhs = DoubleArray(m) { if (indices[it] == j) 1.0 else 0.0 }
        val coeff = LinearAlgebra.solve(matrix, rhs, ctx.backend)
        return ValueFunctional(points, coeff)
    }

    /**
     * Closed-form (explicit) formula of the functional theta_j for an internal index j.
     *
     * Used only as an independent cross-check against the main construction via local
     * biorthogonalization: the two representations must coincide.
     *
     * The functional is based on five points: the nodes x_j, x_{j+3} and the three midpoints of
     * the support intervals. The notation follows the source (docs/REFERENCES.md,
     * section 2): the values of omega_j at these points form the quantities A..E, and the
     * denominator is K1 = C^2 D - B C E - A D E - D E^2.
     *
     * @param j internal index of the functional (boundary j values are not supported: they follow a different branch of the formula)
     * @throws IllegalStateException if the denominator degenerates
     */
    internal fun closedFormInternal(j: Int): ValueFunctional {
        val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2); val xj3 = grid.x(j + 3)
        val a = basis.omega(j, mid(xj, xj1))
        val c = basis.omega(j, mid(xj1, xj2))
        val e = basis.omega(j, mid(xj2, xj3))
        val b = basis.omega(j, xj1)
        val d = basis.omega(j, xj2)
        val k1 = c * c * d - b * c * e - a * d * e - d * e * e
        // The scale of the denominator is the sum of the absolute values of the four terms of K1:
        // the values of omega_j are of order one while their differences are small as powers of h,
        // so an absolute threshold would not be invariant with respect to the grid step.
        // See the KDoc of splines.DEGENERACY_RELATIVE_EPS.
        val k1Scale = cancellationScale(c * c * d, b * c * e, a * d * e, d * e * e)
        check(isSignificant(k1, k1Scale)) {
            "closedFormInternal(j=$j): degenerate denominator K1=$k1, scale=$k1Scale " +
                "(significance lost: threshold $DEGENERACY_RELATIVE_EPS)"
        }
        return ValueFunctional(
            doubleArrayOf(xj, mid(xj, xj1), mid(xj1, xj2), mid(xj2, xj3), xj3),
            doubleArrayOf(e * e / k1, -d * e / k1, (c * d - b * e) / k1, -d * e / k1, e * e / k1),
        )
    }
}
// ----------------------------------------------------------------------------
// xi — de Boor–Fix functionals (value and derivatives)
// Source: Kulikov, Makarov, On de Boor–Fix Type Functionals for Minimal Splines
// (Topics in Classical and Modern Analysis, Springer, 2019, p. 211–225).
// See docs/REFERENCES.md, section 2.
// ----------------------------------------------------------------------------

/**
 * Functional of the form xi(u) = u(node) + cD u'(node) (de Boor–Fix, r = 1, 2).
 * For cD = 0 it reduces to the value u(node) — the boundary functionals u(x_0), u(x_n).
 *
 * @property node node of the functional.
 * @property cD coefficient at the derivative; it has the dimension of length and the order of the grid step.
 */
public class DerivFunctional(public val node: Double, public val cD: Double) : ApproxFunctional {
    override fun apply(f: (Double) -> Double, fD: (Double) -> Double, fDD: (Double) -> Double): Double =
        f(node) + cD * fD(node)

    /**
     * Sum of the absolute values of the representation coefficients `1 + |cD|`.
     *
     * The quantity does not bound the amplification of a data perturbation: the coefficient `cD`
     * is of the order of the step `h`, so `1 + |cD| -> 1` as `h -> 0`, whereas a perturbation eps
     * passed through numerical differentiation is amplified as `eps/h`. The quantity is suitable
     * for diagnosing the representation (the order of the coefficient at the derivative), but not
     * for bounding the effect of data noise.
     */
    override fun absSum(): Double = 1.0 + abs(cD)
}

/**
 * Functional of the form xi^<0>(u) = u(node) + c1 u'(node) + c2 u''(node) (de Boor–Fix, r = 0):
 * it uses the value and the first and second derivatives at a single node.
 *
 * @property node node of the functional.
 * @property c1 coefficient at the first derivative (of order h).
 * @property c2 coefficient at the second derivative (of order h^2).
 */
public class SecondDerivFunctional(public val node: Double, public val c1: Double, public val c2: Double) : ApproxFunctional {
    override fun apply(f: (Double) -> Double, fD: (Double) -> Double, fDD: (Double) -> Double): Double =
        f(node) + c1 * fD(node) + c2 * fDD(node)

    /**
     * Sum of the absolute values of the representation coefficients `1 + |c1| + |c2|`.
     *
     * The quantity does not bound the amplification of a data perturbation for the same reason as
     * in [DerivFunctional]: `c1` is of order `h` and `c2` of order `h^2`, and the sum tends to 1 as
     * the grid is refined, whereas the amplification of a perturbation by the first and second
     * derivatives grows as `1/h` and `1/h^2`. Suitable only for diagnosing the representation.
     */
    override fun absSum(): Double = 1.0 + abs(c1) + abs(c2)
}

/**
 * Family of de Boor–Fix functionals xi_j^{<r>}, r in {0, 1, 2}.
 *
 * All three families are projectors: biorthogonality xi_i(omega_j) = delta_ij holds.
 *
 *  - xi^<1>(u) = u(x_{j+1}) + C1_j u'(x_{j+1}),
 *      C1_j = ((sigma_{j+2}-sigma_{j+1})rho'_{j+2} - (rho_{j+2}-rho_{j+1})sigma'_{j+2}) / W_j;
 *  - xi^<2>(u) = u(x_{j+2}) + C2_j u'(x_{j+2}),
 *      C2_j = ((sigma_{j+2}-sigma_{j+1})rho'_{j+1} - (rho_{j+2}-rho_{j+1})sigma'_{j+1}) / W_j;
 *      where W_j = rho'_{j+2}sigma'_{j+1} - rho'_{j+1}sigma'_{j+2} is the Wronskian;
 *  - xi^<0>(u) = u(x_j) + (N1_j/Delta_j) u'(x_j) + (N2_j/Delta_j) u''(x_j)
 *      — uses the second derivative of its argument, that is, requires u in C^2.
 *
 * The coefficients are computed in the local coordinates of the interval (x_{j+1}, x_{j+2}) (see [splines.LocalFrame]):
 * the formulas for rho, sigma apply to the components psi_1, psi_2, since the first row of the
 * matrix T equals (1, 0, 0) and psi spans the same space. This rules out a loss of significance for
 * a short interval or an interval far from zero; the result does not depend on the choice of
 * coordinates, since the functional is defined by biorthogonality to the basis.
 *
 * The boundary functionals j = -2 and j = n-1 are the values u(x_0) and u(x_n); with this choice
 * biorthogonality holds for all r and all generating systems.
 *
 * @property r order of the functional, one of {0, 1, 2}; r = 1 by default.
 * @param ctx numerical computation context (the family does not use linear algebra).
 * @throws IllegalArgumentException if r is outside {0, 1, 2}, or if on some interval the Wronskian
 *   or the denominator Delta_j is degenerate (see [splines.DEGENERACY_RELATIVE_EPS]).
 */
public class DeBoorFixFunctionals(
    basis: MinimalSplineBasis,
    public val r: Int = 1,
    ctx: NumericsContext = NumericsContext.default(),
) : FunctionalFamily(basis, if (r == 1) "xi" else "xi<$r>", ctx) {
    init { require(r in 0..2) { "DeBoorFix: parameter r must be 0, 1 or 2, got $r" } }
    override val isProjector: Boolean = true
    override val usesDerivative: Boolean = true
    override val usesSecondDerivative: Boolean = (r == 0)
    private val funcs: Array<ApproxFunctional> = Array(n + 2) { buildXi(it - 2) }
    override fun chi(j: Int): ApproxFunctional = funcs[j + 2]

    private fun buildXi(j: Int): ApproxFunctional {
        if (j == -2) return DerivFunctional(grid.x(0), 0.0)
        if (j == n - 1) return DerivFunctional(grid.x(n), 0.0)
        return when (r) {
            0 -> buildXi0(j)
            2 -> buildXi2(j)
            else -> buildXi1(j)
        }
    }

    /** xi^<1>_j: node x_{j+1}, coefficient C1_j. */
    private fun buildXi1(j: Int): ApproxFunctional {
        val x1 = grid.x(j + 1); val x2 = grid.x(j + 2)
        // The components psi_1, psi_2 of the local system of the interval (x_{j+1}, x_{j+2}) and
        // their derivatives: the formula is the same as for rho, sigma, since psi = T phi with the
        // first row of T equal to (1, 0, 0) spans the same space.
        val fr = basis.frame(j + 1)
        val p1 = fr.psi(x1); val p2 = fr.psi(x2)
        val d1 = fr.psiD(x1); val d2 = fr.psiD(x2)
        val rho1 = p1[1]; val rho2 = p2[1]
        val sig1 = p1[2]; val sig2 = p2[2]
        val rhoD1 = d1[1]; val rhoD2 = d2[1]
        val sigD1 = d1[2]; val sigD2 = d2[2]
        // The Wronskian W_j is a difference of two products; its scale is the sum of their absolute
        // values. This makes the check independent of the scale of rho', sigma' themselves (for the
        // system H they grow like cosh, for T they are bounded by one, and on a fine grid the
        // difference is small).
        val denom = rhoD2 * sigD1 - rhoD1 * sigD2
        val denomScale = cancellationScale(rhoD2 * sigD1, rhoD1 * sigD2)
        require(isSignificant(denom, denomScale)) {
            "buildXi1(j=$j): degenerate Wronskian rhoD2*sigD1 - rhoD1*sigD2=$denom, " +
                "scale=$denomScale (significance lost: threshold $DEGENERACY_RELATIVE_EPS)"
        }
        val cD = ((sig2 - sig1) * rhoD2 - (rho2 - rho1) * sigD2) / denom
        return DerivFunctional(x1, cD)
    }

    /** xi^<2>_j: node x_{j+2}, coefficient C2_j (the same denominator W_j, primes taken at x_{j+1}). */
    private fun buildXi2(j: Int): ApproxFunctional {
        val x1 = grid.x(j + 1); val x2 = grid.x(j + 2)
        // The components psi_1, psi_2 of the local system of the interval (x_{j+1}, x_{j+2}) and
        // their derivatives: the formula is the same as for rho, sigma, since psi = T phi with the
        // first row of T equal to (1, 0, 0) spans the same space.
        val fr = basis.frame(j + 1)
        val p1 = fr.psi(x1); val p2 = fr.psi(x2)
        val d1 = fr.psiD(x1); val d2 = fr.psiD(x2)
        val rho1 = p1[1]; val rho2 = p2[1]
        val sig1 = p1[2]; val sig2 = p2[2]
        val rhoD1 = d1[1]; val rhoD2 = d2[1]
        val sigD1 = d1[2]; val sigD2 = d2[2]
        // The same Wronskian as in buildXi1: the scale is the sum of the absolute values of the two products.
        val denom = rhoD2 * sigD1 - rhoD1 * sigD2
        val denomScale = cancellationScale(rhoD2 * sigD1, rhoD1 * sigD2)
        require(isSignificant(denom, denomScale)) {
            "buildXi2(j=$j): degenerate Wronskian rhoD2*sigD1 - rhoD1*sigD2=$denom, " +
                "scale=$denomScale (significance lost: threshold $DEGENERACY_RELATIVE_EPS)"
        }
        val cD = ((sig2 - sig1) * rhoD1 - (rho2 - rho1) * sigD1) / denom
        return DerivFunctional(x2, cD)
    }

    /**
     * xi^<0>_j: node x_j, coefficients N1_j/Delta_j (at u') and N2_j/Delta_j (at u'') following the
     * formulas of the source; Delta_j = W_j (rho'_j sigma''_j - rho''_j sigma'_j).
     */
    private fun buildXi0(j: Int): ApproxFunctional {
        val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2)
        // The values of psi_1, psi_2 and their derivatives at three nodes in the local coordinates
        // of the interval (x_{j+1}, x_{j+2}); the node x_j is at distance h_j from the origin.
        val fr = basis.frame(j + 1)
        val pj = fr.psi(xj); val dj = fr.psiD(xj); val ddj = fr.psiDD(xj)
        val pj1 = fr.psi(xj1); val dj1 = fr.psiD(xj1)
        val pj2 = fr.psi(xj2); val dj2 = fr.psiD(xj2)
        val rj = pj[1]; val sj = pj[2]
        val rDj = dj[1]; val sDj = dj[2]
        val rDDj = ddj[1]; val sDDj = ddj[2]
        val rj1 = pj1[1]; val sj1 = pj1[2]
        val rDj1 = dj1[1]; val sDj1 = dj1[2]
        val rj2 = pj2[1]; val sj2 = pj2[2]
        val rDj2 = dj2[1]; val sDj2 = dj2[2]

        // Delta_j is the product of two 2x2 minors and degenerates exactly when one of the factors
        // degenerates. The factors are checked separately, each against its own scale: this way the
        // diagnostic points at the cause and the criterion does not depend on the product of the scales.
        val wronskian12 = rDj1 * sDj2 - rDj2 * sDj1
        val wronskian12Scale = cancellationScale(rDj1 * sDj2, rDj2 * sDj1)
        require(isSignificant(wronskian12, wronskian12Scale)) {
            "buildXi0(j=$j): degenerate Wronskian W_j=$wronskian12, scale=$wronskian12Scale " +
                "(significance lost: threshold $DEGENERACY_RELATIVE_EPS)"
        }
        val curvature = rDj * sDDj - rDDj * sDj
        val curvatureScale = cancellationScale(rDj * sDDj, rDDj * sDj)
        require(isSignificant(curvature, curvatureScale)) {
            "buildXi0(j=$j): degenerate denominator rho'sigma''-rho''sigma'=$curvature, scale=$curvatureScale " +
                "(significance lost: threshold $DEGENERACY_RELATIVE_EPS)"
        }
        val delta = wronskian12 * curvature
        // Checking the factors does not control the product itself: two significant but small
        // factors (of order 1e-200) give delta == 0.0 by underflow, two large ones give Inf by
        // overflow; in both cases n1/delta would be non-finite. `require` is used, as for the
        // factors: every case of invalid input ends with an IllegalArgumentException.
        require(delta.isFinite() && delta != 0.0) {
            "buildXi0(j=$j): Delta_j=$delta is unusable as a denominator (underflow/overflow of the product " +
                "of individually significant factors): W_j=$wronskian12, " +
                "rho'sigma''-rho''sigma'=$curvature"
        }
        val n1 = (rj1 * sDj1 - rDj1 * sj1) * (rDDj * sDj2 - rDj2 * sDDj) +
            (rDj1 * sDj2 - rDj2 * sDj1) * (rDDj * sj - rj * sDDj) +
            (rj2 * sDj2 - rDj2 * sj2) * (rDj1 * sDDj - rDDj * sDj1)
        val n2 = (rj1 * sDj1 - rDj1 * sj1) * (rDj2 * sDj - rDj * sDj2) +
            (rDj1 * sDj2 - rDj2 * sDj1) * (rj * sDj - rDj * sj) +
            (rj2 * sDj2 - rDj2 * sj2) * (rDj * sDj1 - rDj1 * sDj)
        // Finiteness of the denominator is not enough: for a subnormal non-zero delta (of order
        // 1e-310) the quotients n1/delta, n2/delta overflow, so the result itself is checked.
        val c1 = n1 / delta
        val c2 = n2 / delta
        require(c1.isFinite() && c2.isFinite()) {
            "buildXi0(j=$j): non-finite coefficients: N1/Delta=$c1, N2/Delta=$c2 " +
                "(N1=$n1, N2=$n2, Delta_j=$delta — overflow when dividing by a subnormal denominator)"
        }
        return SecondDerivFunctional(xj, c1, c2)
    }
}

// ----------------------------------------------------------------------------
// xitilde — discretized de Boor–Fix functionals (without derivatives)
// ----------------------------------------------------------------------------

/**
 * Discretized de Boor–Fix functionals xitilde^{<r>}_j, r in {1, 2}, which do not use the
 * derivative.
 *
 * The derivative f'(x_k) in xi^{<1>}, xi^{<2>} is replaced by the central divided difference of
 * nodal values f'(x_k) ≈ (f(x_{k+1}) - f(x_{k-1})) / (x_{k+1} - x_{k-1}), which is also valid on
 * non-uniform grids:
 *   xitilde^{<1>}_j(f) = f(x_{j+1}) + w1_j (f(x_{j+2}) - f(x_j)) / (x_{j+2} - x_j),
 *   xitilde^{<2>}_j(f) = f(x_{j+2}) + w2_j (f(x_{j+3}) - f(x_{j+1})) / (x_{j+3} - x_{j+1}),
 * where w1_j, w2_j are the coefficients at the derivative from [DeBoorFixFunctionals]
 * ([DerivFunctional.cD]). The functionals are represented as [ValueFunctional] (usesDerivative = false).
 *
 * The operator P_xitilde is a quasi-interpolant, not a projector (isProjector = false): replacing
 * the derivative by a difference destroys exact biorthogonality in general. On span{1, rho, sigma}
 * the error is of order O(h^2 |f''|); with multiple boundary nodes the central difference becomes
 * one-sided at the boundary, and the order drops in the boundary layer. The construction is our own
 * (docs/REFERENCES.md, section 2).
 *
 * The boundary functionals j = -2, n-1 are the values f(x_0), f(x_n), as in theta and xi.
 *
 * The variant xitilde^{<0>} is not implemented: it would rely on f'' at the left end of the support,
 * where the spline and its first derivative vanish, and replacing the second derivative by a
 * difference does not reproduce the generating system even as h → 0.
 *
 * @property r order of the functional, one of {1, 2}; r = 1 by default.
 * @param ctx numerical computation context; passed to the nested [DeBoorFixFunctionals] family.
 * @throws IllegalArgumentException if r is outside {1, 2} or the nested xi family cannot be built.
 */
public class DiscreteDeBoorFixFunctionals(
    basis: MinimalSplineBasis,
    public val r: Int = 1,
    ctx: NumericsContext = NumericsContext.default(),
) : FunctionalFamily(basis, if (r == 1) "xitilde" else "xitilde<$r>", ctx) {
    init { require(r in 1..2) { "DiscreteDeBoorFix: parameter r must be 1 or 2, got $r" } }
    override val isProjector: Boolean = false
    override val usesDerivative: Boolean = false
    // The context is passed to the nested family so that `raw.ctx` matches the context of the wrapper.
    private val raw = DeBoorFixFunctionals(basis, r, ctx)
    private val funcs: Array<ApproxFunctional> = Array(n + 2) { buildXiTilde(it - 2) }
    override fun chi(j: Int): ApproxFunctional = funcs[j + 2]

    private fun buildXiTilde(j: Int): ApproxFunctional {
        if (j == -2) return ValueFunctional(doubleArrayOf(grid.x(0)), doubleArrayOf(1.0))
        if (j == n - 1) return ValueFunctional(doubleArrayOf(grid.x(n)), doubleArrayOf(1.0))
        // The factor at the derivative is the same A^{<r>}_j as in the original xi.
        val w = (raw.chi(j) as DerivFunctional).cD
        // The support node and its neighbours for the central difference around it.
        val (node, left, right) = when (r) {
            2 -> Triple(grid.x(j + 2), grid.x(j + 1), grid.x(j + 3))
            else -> Triple(grid.x(j + 1), grid.x(j), grid.x(j + 2))
        }
        // The step of the divided difference is the difference of node coordinates; its scale is
        // |right| + |left|. The difference itself is small as h, so an absolute threshold would not
        // tell a fine grid from coincident nodes on an interval far from zero.
        // See the KDoc of splines.DEGENERACY_RELATIVE_EPS.
        val denom = right - left
        val denomScale = cancellationScale(right, left)
        require(isSignificant(denom, denomScale)) {
            "buildXiTilde(j=$j,r=$r): degenerate divided-difference step x=$right - x=$left = $denom, " +
                "scale=$denomScale (significance lost: threshold $DEGENERACY_RELATIVE_EPS)"
        }
        // f(node) + w (f(right) - f(left))/denom as a combination of values.
        return ValueFunctional(
            doubleArrayOf(left, node, right),
            doubleArrayOf(-w / denom, 1.0, w / denom),
        )
    }
}

// ----------------------------------------------------------------------------
// mu — averaging functionals
// Source: Kulikov, Makarov, Construction of Approximation Functionals for
// Minimal Splines (Journal of Mathematical Sciences, 2022, vol. 262, no. 1,
// p. 84–98). See docs/REFERENCES.md, section 2.
// ----------------------------------------------------------------------------

/**
 * Averaging functionals mu_j(f) = a_j f(y_{j-1}) + b_j f(y_j) + c_j f(y_{j+1}).
 *
 * Nodes of the auxiliary grid: y_j = x_{j+1} + theta (x_{j+2} - x_{j+1}), with theta = 1/2 by
 * default (the interval midpoints). The coefficients are determined from the condition of exactness
 * on span{1, rho, sigma}: the system
 *   [1, 1, 1; rho(y_{j-1}), rho(y_j), rho(y_{j+1}); sigma(...)] (a, b, c)^T = a^N_j
 * is solved, where a^N_j is the same approximation-relation vector that the basis itself builds.
 *
 * The operator P_mu is a quasi-interpolant, not a projector (P_mu^2 != P_mu); it is exact on
 * span{1, rho, sigma}.
 *
 * Reference special case: for the polynomial system B with theta = 1/2 on a uniform grid the
 * formula reduces to -1/8 (f(y_{j-1}) - 10 f(y_j) + f(y_{j+1})).
 *
 * @property theta placement parameter of the auxiliary grid nodes, in (0, 1).
 * @param ctx numerical computation context for the 3×3 linear systems.
 */
public class AveragingFunctionals(
    basis: MinimalSplineBasis,
    public val theta: Double = 0.5,
    ctx: NumericsContext = NumericsContext.default(),
) : FunctionalFamily(basis, "mu", ctx) {
    override val isProjector: Boolean = false
    override val usesDerivative: Boolean = false
    private val funcs: Array<ApproxFunctional> = Array(n + 2) { buildMu(it - 2) }
    override fun chi(j: Int): ApproxFunctional = funcs[j + 2]

    /** y_j per (net_Y): boundary y_{-2}=x_0, y_{n-1}=x_n; internal — x_{j+1}+theta(x_{j+2}-x_{j+1}). */
    private fun yNode(j: Int): Double = when (j) {
        -2 -> grid.x(0)
        n - 1 -> grid.x(n)
        else -> grid.x(j + 1) + theta * (grid.x(j + 2) - grid.x(j + 1))
    }

    private fun buildMu(j: Int): ApproxFunctional {
        if (j == -2) return ValueFunctional(doubleArrayOf(grid.x(0)), doubleArrayOf(1.0))
        if (j == n - 1) return ValueFunctional(doubleArrayOf(grid.x(n)), doubleArrayOf(1.0))
        val ym = yNode(j - 1); val y0 = yNode(j); val yp = yNode(j + 1)
        val ys = doubleArrayOf(ym, y0, yp)
        // The system (phi(y_{j-1}) | phi(y_j) | phi(y_{j+1})) mu = a_j is written in the local
        // coordinates of the interval (x_{j+1}, x_{j+2}) — the middle interval of the support of
        // omega_j: multiplying both sides on the left by T_{j+1} does not change the solution, while
        // the conditioning of the matrix stops depending on the grid step and the position of the
        // interval. All three points y_q lie within two steps of x_{j+1}, where psi = O(1). The
        // right-hand side T_{j+1} a_j is obtained by the same approximation-relation formula as
        // the columns T_k M_k in the basis.
        val frame = basis.frame(j + 1)
        val cols = Array(3) { q -> frame.psi(ys[q]) }
        val matrix = DenseMatrix.build(3, 3) { r, c -> cols[c][r] }
        val coeff = LinearAlgebra.solve(matrix, basis.computeA(j, frame.psi, frame.psiD), ctx.backend)
        return ValueFunctional(ys, coeff)
    }
}

// ----------------------------------------------------------------------------
// lambda — three-point functionals
// Source: Kulikov, Makarov (Journal of Mathematical Sciences, 2022, vol. 262,
// no. 1, p. 84–98). See docs/REFERENCES.md, section 2.
// ----------------------------------------------------------------------------

/**
 * Three-point functionals lambda_j(f) based on the points x_{j+1}, x_{j+3/2}, x_{j+2},
 * where x_{j+3/2} = x_{j+1} + thetaHat (x_{j+2} - x_{j+1}).
 *
 * They are implemented via a local approximation on the interval I = [x_{j+1}, x_{j+2}]: exactly
 * three splines omega_{j-1}, omega_j, omega_{j+1} are active on it; the system at the three points
 * is solved and lambda_j(f) is taken as the coefficient at omega_j.
 *
 * The operator P_lambda is a quasi-interpolant, not a projector (P_lambda^2 != P_lambda); it is
 * exact on span{1, rho, sigma}.
 *
 * Reference special case: for the system B with thetaHat = 1/2 the formula reduces
 * to -1/2 (f(x_{j+1}) - 4 f(x_{j+3/2}) + f(x_{j+2})).
 *
 * @property thetaHat placement parameter of the middle point, in (0, 1); 1/2 by default.
 * @param ctx numerical computation context for the 3×3 linear systems.
 */
public class ThreePointFunctionals(
    basis: MinimalSplineBasis,
    public val thetaHat: Double = 0.5,
    ctx: NumericsContext = NumericsContext.default(),
) : FunctionalFamily(basis, "lambda", ctx) {
    override val isProjector: Boolean = false
    override val usesDerivative: Boolean = false
    private val funcs: Array<ApproxFunctional> = Array(n + 2) { buildLambda(it - 2) }
    override fun chi(j: Int): ApproxFunctional = funcs[j + 2]

    private fun buildLambda(j: Int): ApproxFunctional {
        if (j == -2) return ValueFunctional(doubleArrayOf(grid.x(0)), doubleArrayOf(1.0))
        if (j == n - 1) return ValueFunctional(doubleArrayOf(grid.x(n)), doubleArrayOf(1.0))
        val x1 = grid.x(j + 1); val x2 = grid.x(j + 2)
        val xMid = x1 + thetaHat * (x2 - x1)
        val points = doubleArrayOf(x1, xMid, x2)
        val active = intArrayOf(j - 1, j, j + 1) // active on (x_{j+1},x_{j+2})
        // M[p][slot] = omega_active[slot](point_p); from M c = fvals it follows that lambda_j = c[1].
        // The coefficients coeff_p = (M^{-1})[1][p] solve M^T coeff = e_1 (the second row of the inverse).
        val mTrans = DenseMatrix.build(3, 3) { i, p -> basis.omega(active[i], points[p]) }
        val e1 = doubleArrayOf(0.0, 1.0, 0.0)
        val coeff = LinearAlgebra.solve(mTrans, e1, ctx.backend)
        return ValueFunctional(points, coeff)
    }
}
