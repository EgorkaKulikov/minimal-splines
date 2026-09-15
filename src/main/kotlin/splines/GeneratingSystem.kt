package splines

import kotlin.math.cos
import kotlin.math.cosh
import kotlin.math.min
import kotlin.math.sin
import kotlin.math.sinh

// ============================================================================
// Generating vector function phi(t) = (1, rho(t), sigma(t))^T
// ============================================================================

/**
 * Local representation of the generating system on a grid interval: psi(t) = T phi(t), where T is a
 * non-singular 3×3 matrix chosen so that the components of psi and of its derivatives are of order one
 * on the interval. After left multiplication by T the approximation relation
 * sum_j a_j omega_j(t) = phi(t) turns into the equivalent (T M_k) omega(t) = psi(t) with the same
 * solution omega(t), so the basis functions do not depend on the choice of representation, while the
 * conditioning of the matrix T M_k stops depending on the number of intervals and on the position and
 * scale of the interval.
 *
 * The values [psi], [psiD], [psiDD] are computed directly in the local variable rather than as the
 * product of T and phi(t): the latter would produce mutual cancellation of large terms.
 *
 * @property psi local vector function psi(t) = T phi(t);
 * @property psiD its derivative psi'(t) = T phi'(t);
 * @property psiDD the second derivative psi''(t) = T phi''(t);
 * @property det determinant of the matrix T.
 */
public class LocalFrame(
    public val psi: (Double) -> DoubleArray,
    public val psiD: (Double) -> DoubleArray,
    public val psiDD: (Double) -> DoubleArray,
    public val det: Double,
)

/**
 * Generating vector function phi(t) = (1, rho(t), sigma(t))^T and its derivatives up to the second
 * order. The component phi_0 = 1 makes the basis a partition of unity.
 *
 * The optional parameter localFrameFactory defines the local representation of the system on an
 * interval (see [localFrame]); it is defined for the built-in systems [B], [H], [T]. A system without
 * it is used in global coordinates.
 *
 * @property name short name of the system (B, H, T or a user-defined one).
 * @property rho second component rho(t).
 * @property sigma third component sigma(t).
 * @property rhoD derivative rho'(t).
 * @property sigmaD derivative sigma'(t).
 * @property rhoDD second derivative rho''(t).
 * @property sigmaDD second derivative sigma''(t).
 */
public class GeneratingSystem(
    public val name: String,
    public val rho: (Double) -> Double,
    public val sigma: (Double) -> Double,
    public val rhoD: (Double) -> Double,
    public val sigmaD: (Double) -> Double,
    public val rhoDD: (Double) -> Double,
    public val sigmaDD: (Double) -> Double,
    private val localFrameFactory: ((c: Double, h: Double) -> LocalFrame)? = null,
) {
    /** phi(t) = (1, rho(t), sigma(t)). */
    public fun phi(t: Double): DoubleArray = doubleArrayOf(1.0, rho(t), sigma(t))

    /** phi'(t) = (0, rho'(t), sigma'(t)). */
    public fun phiD(t: Double): DoubleArray = doubleArrayOf(0.0, rhoD(t), sigmaD(t))

    /** phi''(t) = (0, rho''(t), sigma''(t)). */
    public fun phiDD(t: Double): DoubleArray = doubleArrayOf(0.0, rhoDD(t), sigmaDD(t))

    /** Wronskian det(phi(t), phi'(t), phi''(t)); non-zero for a non-degenerate system. */
    public fun wronskian(t: Double): Double {
        val a = phi(t); val b = phiD(t); val c = phiDD(t)
        return a[0] * (b[1] * c[2] - b[2] * c[1]) +
            a[1] * (b[2] * c[0] - b[0] * c[2]) +
            a[2] * (b[0] * c[1] - b[1] * c[0])
    }

    /**
     * Local representation of the system on the interval [c, c + h] (see [LocalFrame]).
     *
     * For systems closed under a shift of the argument (solutions of linear homogeneous differential
     * equations with constant coefficients, including [B], [H], [T]), the representation is built by
     * shifting to the point c and scaling by the step h. For a system that is not closed under a shift
     * of the argument the local representation is unavailable: the global representation with T = I is
     * used, that is, psi = phi.
     */
    public fun localFrame(c: Double, h: Double): LocalFrame =
        localFrameFactory?.invoke(c, h) ?: LocalFrame(::phi, ::phiD, ::phiDD, 1.0)

    /** Built-in generating systems. */
    public companion object {
        /** Polynomial phi^B(t) = (1, t, t^2)^T. */
        public val B: GeneratingSystem = GeneratingSystem(
            name = "B",
            rho = { t -> t }, sigma = { t -> t * t },
            rhoD = { 1.0 }, sigmaD = { t -> 2.0 * t },
            rhoDD = { 0.0 }, sigmaDD = { 2.0 },
            localFrameFactory = ::polynomialFrame,
        )

        /** Hyperbolic phi^H(t) = (1, sinh t, cosh t)^T. */
        public val H: GeneratingSystem = GeneratingSystem(
            name = "H",
            rho = { t -> Math.sinh(t) }, sigma = { t -> Math.cosh(t) },
            rhoD = { t -> Math.cosh(t) }, sigmaD = { t -> Math.sinh(t) },
            rhoDD = { t -> Math.sinh(t) }, sigmaDD = { t -> Math.cosh(t) },
            localFrameFactory = ::hyperbolicFrame,
        )

        /** Trigonometric phi^T(t) = (1, sin t, cos t)^T. */
        public val T: GeneratingSystem = GeneratingSystem(
            name = "T",
            rho = { t -> Math.sin(t) }, sigma = { t -> Math.cos(t) },
            rhoD = { t -> Math.cos(t) }, sigmaD = { t -> -Math.sin(t) },
            rhoDD = { t -> -Math.sin(t) }, sigmaDD = { t -> -Math.cos(t) },
            localFrameFactory = ::trigonometricFrame,
        )
    }
}

/**
 * Local representation of the system B: psi(t) = (1, s, s^2), s = (t - c)/h.
 * The matrix T has rows (1, 0, 0), (-c/h, 1/h, 0), (c^2/h^2, -2c/h^2, 1/h^2); det T = 1/h^3.
 */
private fun polynomialFrame(c: Double, h: Double): LocalFrame {
    val ih = 1.0 / h
    return LocalFrame(
        psi = { t -> val s = (t - c) * ih; doubleArrayOf(1.0, s, s * s) },
        psiD = { t -> doubleArrayOf(0.0, ih, 2.0 * (t - c) * ih * ih) },
        psiDD = { doubleArrayOf(0.0, 0.0, 2.0 * ih * ih) },
        det = ih * ih * ih,
    )
}

/**
 * Scale of the local representation for the systems H and T: the step h, but not greater than one —
 * the characteristic scale of these systems; for h > 1 the components sinh(t - c), 1 - cos(t - c) are
 * already of order one.
 */
private fun boundedScale(h: Double): Double = min(h, 1.0)

/**
 * Local representation of the system H: psi(t) = (1, sinh(u)/l, (cosh(u) - 1)/l^2), u = t - c,
 * l = min(h, 1). The third component is computed as 2 sinh^2(u/2)/l^2 without loss of significance.
 * The matrix T has rows (1, 0, 0), (0, cosh c/l, -sinh c/l), (-1/l^2, -sinh c/l^2, cosh c/l^2);
 * det T = 1/l^3.
 */
private fun hyperbolicFrame(c: Double, h: Double): LocalFrame {
    val il = 1.0 / boundedScale(h)
    val il2 = il * il
    return LocalFrame(
        psi = { t -> val u = t - c; val sh = sinh(0.5 * u); doubleArrayOf(1.0, sinh(u) * il, 2.0 * sh * sh * il2) },
        psiD = { t -> val u = t - c; doubleArrayOf(0.0, cosh(u) * il, sinh(u) * il2) },
        psiDD = { t -> val u = t - c; doubleArrayOf(0.0, sinh(u) * il, cosh(u) * il2) },
        det = il * il2,
    )
}

/**
 * Local representation of the system T: psi(t) = (1, sin(u)/l, (1 - cos(u))/l^2), u = t - c,
 * l = min(h, 1). The third component is computed as 2 sin^2(u/2)/l^2 without loss of significance.
 * The matrix T has rows (1, 0, 0), (0, cos c/l, -sin c/l), (1/l^2, -sin c/l^2, -cos c/l^2);
 * det T = -1/l^3.
 */
private fun trigonometricFrame(c: Double, h: Double): LocalFrame {
    val il = 1.0 / boundedScale(h)
    val il2 = il * il
    return LocalFrame(
        psi = { t -> val u = t - c; val s = sin(0.5 * u); doubleArrayOf(1.0, sin(u) * il, 2.0 * s * s * il2) },
        psiD = { t -> val u = t - c; doubleArrayOf(0.0, cos(u) * il, sin(u) * il2) },
        psiDD = { t -> val u = t - c; doubleArrayOf(0.0, -sin(u) * il, cos(u) * il2) },
        det = -il * il2,
    )
}
