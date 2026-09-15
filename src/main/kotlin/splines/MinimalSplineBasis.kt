package splines

import numerics.Conditioning
import numerics.DenseMatrix
import numerics.NumericsContext

// ============================================================================
// Minimal spline basis: omega(t) = M_k^{-1} phi(t) in the local coordinates of the interval
// ============================================================================

/**
 * Basis of quadratic minimal splines {omega_j}_{j=-2}^{n-1} on a grid with triple boundary nodes.
 *
 * On the interval (x_k, x_{k+1}) the values are defined by the approximation relation
 * sum_j a_j omega_j(t) = phi(t), that is, omega(t) = M_k^{-1} phi(t), where M_k = (a_{k-2}|a_{k-1}|a_k).
 * The computation is carried out in the local representation of the generating system psi_k = T_k phi
 * (see [GeneratingSystem.localFrame]): left multiplication by the non-singular matrix T_k yields the
 * equivalent system (T_k M_k) omega(t) = psi_k(t) with the same solution, so the basis functions do
 * not depend on the choice of representation. The columns T_k a_j are computed by the formula
 * [computeA] written for psi_k instead of phi: the vector a_j is covariant under the substitution
 * phi -> T phi (a_j^{T phi} = T a_j^{phi}), since the coefficient at phi'(x_{j+1}) is the ratio of two
 * dot products with one and the same normal vector and is unchanged by a linear substitution.
 * The components of psi_k are of order one on the interval, and the condition number of T_k M_k
 * depends neither on the number of intervals nor on the position and scale of the interval.
 *
 * The inverse matrices are computed by the BLAS/LAPACK implementation taken from the context ctx; the
 * reliability of the inversion is controlled by a condition number estimate for T_k M_k
 * (see [MAX_CONDITION]).
 *
 * @property sys generating system phi = (1, rho, sigma).
 * @property grid grid with multiple boundary nodes.
 * @param ctx numerics context (BLAS/LAPACK implementation).
 * @throws IllegalArgumentException if on some interval the generating system overflows in double, or
 *   the approximation-relation matrix is singular, or its condition number exceeds [MAX_CONDITION].
 */
public class MinimalSplineBasis(public val sys: GeneratingSystem, public val grid: Grid, ctx: NumericsContext = NumericsContext.default()) {
    /** Constants of the degeneracy criterion. */
    public companion object {
        /**
         * Largest admissible condition number of the approximation-relation matrix in the local
         * coordinates of the interval; for a larger value fewer than eight significant digits survive in
         * the inverse matrix, and the basis is considered degenerate on that interval.
         */
        public const val MAX_CONDITION: Double = 1e8
    }

    /** Number of grid intervals; the basis consists of n + 2 functions omega_{-2}, ..., omega_{n-1}. */
    public val n: Int = grid.n

    /** Local representations of the generating system on the intervals (x_k, x_{k+1}), k = 0..n-1. */
    private val frames: Array<LocalFrame> = Array(n) { k -> sys.localFrame(grid.x(k), grid.x(k + 1) - grid.x(k)) }

    /** Local representation psi_k = T_k phi for the interval (x_k, x_{k+1}); used by the functional families. */
    internal fun frame(k: Int): LocalFrame = frames[k]

    /** Inverse matrices (T_k M_k)^{-1} stored column-wise: the entry (slot, p) is kept in data[slot + 3 p]. */
    private val invM: Array<DoubleArray> = Array(n) { k -> invertApproximationMatrix(k, ctx) }

    /** Approximation-relation matrix M_k = (a_{k-2}|a_{k-1}|a_k) in the global coordinates of phi. */
    internal fun approximationMatrix(k: Int): DenseMatrix {
        val cols = Array(3) { j -> computeA(k - 2 + j) }
        return DenseMatrix.build(3, 3) { i, j -> cols[j][i] }
    }

    /** Matrix T_k M_k = (T_k a_{k-2}|T_k a_{k-1}|T_k a_k) in the local coordinates of interval k; it is the one that gets inverted. */
    internal fun localApproximationMatrix(k: Int): DenseMatrix {
        val frame = frames[k]
        requireFiniteFrame(k, frame)
        val cols = Array(3) { j -> computeA(k - 2 + j, frame.psi, frame.psiD) }
        return DenseMatrix.build(3, 3) { i, j -> cols[j][i] }
    }

    /**
     * Checks that the generating system is representable in double on interval k: the columns a_{k-2},
     * a_{k-1}, a_k are built from the values of psi_k, psi_k' at the nodes x_{k-1}, ..., x_{k+2}, and all
     * of them must be finite. For the systems H and T with the local scale l = min(h, 1) the values
     * sinh(u/l) already overflow for u/l > 710, that is, for a step h > 355 in the case of H; without
     * this check the overflow would show up as a NaN in the message about a singular approximation
     * relation.
     */
    private fun requireFiniteFrame(k: Int, frame: LocalFrame) {
        for (m in k - 1..k + 2) {
            val x = grid.x(m)
            require(allFinite(frame.psi(x)) && allFinite(frame.psiD(x))) {
                "Generating system ${sys.name} overflows on interval $k = [${grid.x(k)}, ${grid.x(k + 1)}]: " +
                    "the values of psi_k or psi_k' at the node x_$m = $x are non-finite; reduce the grid step or the length of the interval"
            }
        }
    }

    private fun invertApproximationMatrix(k: Int, ctx: NumericsContext): DoubleArray {
        val m = localApproximationMatrix(k)
        val cond = Conditioning.conditionEstimate(m, ctx).valueOrNull()
            ?: throw IllegalArgumentException(
                "Approximation relation matrix on interval $k is numerically singular: " +
                    "the condition number estimate is unreliable",
            )
        if (cond > MAX_CONDITION) {
            throw IllegalArgumentException(
                "Approximation relation matrix on interval $k is ill-conditioned: " +
                    "condition number $cond exceeds MAX_CONDITION = $MAX_CONDITION",
            )
        }
        val inv = Conditioning.inverse(m, ctx.backend)
            ?: throw IllegalArgumentException(
                "Approximation relation matrix on interval $k is numerically singular",
            )
        return inv.data
    }

    /**
     * Vector a_j of the approximation relation on (x_{j+1}, x_{j+2}), j = -2..n-1, in the global
     * coordinates of phi.
     *
     * `internal` rather than `private`: the same vector a^N_j is needed by the family of averaging
     * functionals mu (`AveragingFunctionals`). It takes no part in building the basis: the columns of
     * T_k M_k are computed through the local representation [frame]. For the systems H and T on
     * intervals with |t| > 700 the global values of sinh, cosh overflow, and the method is inapplicable
     * (see [requireFiniteFrame]).
     */
    internal fun computeA(j: Int): DoubleArray = computeA(j, sys::phi, sys::phiD)

    /**
     * Vector a_j for the generating vector function [phi] with derivative [phiD]: the direction vector
     * of the intersection of the planes span{phi(x_{j+1}), phi'(x_{j+1})} and
     * span{phi(x_{j+2}), phi'(x_{j+2})}, normalised by the condition
     * a_j = phi(x_{j+1}) - coef phi'(x_{j+1}). The coefficient coef is the ratio of the dot products
     * with the common normal phi(x_{j+2}) × phi'(x_{j+2}), so under the substitution phi -> T phi the
     * vector turns into T a_j; this allows the columns of T_k M_k to be computed by the same formula.
     */
    internal fun computeA(j: Int, phi: (Double) -> DoubleArray, phiD: (Double) -> DoubleArray): DoubleArray {
        val xj1 = grid.x(j + 1)
        val phiJ1 = phi(xj1)
        if (grid.isCoincident(j + 1)) return phiJ1 // triple node at the boundary: x_{j+1} = x_{j+2}
        val xj2 = grid.x(j + 2)
        val phiDJ1 = phiD(xj1)
        val dJ2 = cross(phi(xj2), phiD(xj2))
        // The denominator is a dot product; its scale is given by the sum of the absolute values of the
        // componentwise products (the magnitude before mutual cancellation). An absolute threshold is
        // not applicable here: denom ~ h, that is, it depends on the scale of the interval.
        val denom = dot(dJ2, phiDJ1)
        val denomScale = dot3Scale(dJ2, phiDJ1)
        require(denomScale.isFinite()) {
            "computeA(j=$j): products of the values of phi and phi' at the nodes x_{j+1}=$xj1, x_{j+2}=$xj2 " +
                "overflow (scale $denomScale); reduce the grid step or the length of the interval"
        }
        require(isSignificant(denom, denomScale)) {
            "computeA(j=$j): degenerate approximation relation, dot(dJ2, phiDJ1)=$denom, " +
                "scale=$denomScale (significance lost: threshold $DEGENERACY_RELATIVE_EPS)"
        }
        val coef = dot(dJ2, phiJ1) / denom
        return doubleArrayOf(
            phiJ1[0] - coef * phiDJ1[0],
            phiJ1[1] - coef * phiDJ1[1],
            phiJ1[2] - coef * phiDJ1[2],
        )
    }
    /**
     * Index of the grid interval k with x_k <= t < x_{k+1} (for t = b the value n-1 is returned).
     *
     * The single place where membership of the point in the interval is checked: [interval],
     * [evalSpline], [evalSplineDeriv], [evalSplineDeriv2] all go through it. The omega* methods with an
     * outside point do not reach it: they are preceded by a cut-off against the support
     * [x_j, x_{j+3}], which lies inside the interval.
     *
     * @throws IllegalArgumentException if t lies strictly outside [Grid.a], [Grid.b].
     */
    private fun intervalOf(t: Double): Int {
        // The condition is written with negations rather than as `t in grid.a..grid.b`: for t = NaN both
        // comparisons are false, so NaN passes through and propagates into the result (phi(NaN) = NaN);
        // the range form would have turned NaN into an exception.
        require(!(t < grid.a) && !(t > grid.b)) {
            val side = if (t < grid.a) "to the left of a by ${grid.a - t}" else "to the right of b by ${t - grid.b}"
            "MinimalSplineBasis: point t=$t is outside the grid interval [${grid.a}, ${grid.b}] ($side). " +
                "The spline is defined only on the grid interval, and outside it the value is not extrapolated. " +
                "If the point comes from a finite-difference stencil or a quadrature with t > b, " +
                "the argument must be clamped to the interval."
        }
        // Binary search (the breakpoints x_0..x_n are increasing): the largest k in [0,n-1]
        // with x_k <= t, clamped at the ends. The result matches a linear search:
        // t < x_1 -> 0; t >= x_{n-1} -> n-1; otherwise x_lo <= t < x_{lo+1}.
        if (t < grid.x(1)) return 0
        if (t >= grid.x(n - 1)) return n - 1
        var lo = 1
        var hi = n - 1
        while (hi - lo > 1) {
            val mid = (lo + hi) / 2
            if (grid.x(mid) <= t) lo = mid else hi = mid
        }
        return lo
    }

    /**
     * Index of the grid interval containing t (public access).
     *
     * @throws IllegalArgumentException if t is strictly outside the grid interval (see [intervalOf]).
     */
    public fun interval(t: Double): Int = intervalOf(t)

    /** Checks the index j of a basis spline, which must lie in [-2, n-1]. */
    private fun requireIndex(j: Int) {
        require(j in -2..n - 1) { "Spline index j must lie in [-2, ${n - 1}], got $j" }
    }

    /** Checks the length of the coefficient vector: one per basis spline omega_{-2}, ..., omega_{n-1}. */
    private fun requireCoefficients(c: DoubleArray) {
        require(c.size == n + 2) { "Coefficient vector must have length n + 2 = ${n + 2}, got ${c.size}" }
    }

    /**
     * The three active values omega_{k-2},omega_{k-1},omega_k at the point t (one (T_k M_k)^{-1} psi_k(t)).
     *
     * @throws IllegalArgumentException if k is outside [0, n-1].
     */
    public fun activeOmega(k: Int, t: Double): DoubleArray {
        require(k in 0..n - 1) { "Interval index k must lie in [0, ${n - 1}], got $k" }
        val inv = invM[k]
        val p = frames[k].psi(t)
        return doubleArrayOf(
            inv[0] * p[0] + inv[3] * p[1] + inv[6] * p[2],
            inv[1] * p[0] + inv[4] * p[1] + inv[7] * p[2],
            inv[2] * p[0] + inv[5] * p[1] + inv[8] * p[2],
        )
    }

    /**
     * Value of omega_j(t), j in [-2, n-1]. Support [x_j, x_{j+3}]; outside it the value is 0.
     *
     * @throws IllegalArgumentException if j is outside [-2, n-1].
     */
    public fun omega(j: Int, t: Double): Double {
        requireIndex(j)
        if (t < grid.x(j) || t > grid.x(j + 3)) return 0.0
        val k = intervalOf(t)
        val slot = j - (k - 2)
        if (slot < 0 || slot > 2) return 0.0
        val inv = invM[k]
        val p = frames[k].psi(t)
        return inv[slot] * p[0] + inv[slot + 3] * p[1] + inv[slot + 6] * p[2]
    }

    /**
     * Derivative omega_j'(t) (psi_k is replaced by psi_k'). Needed by the xi functionals.
     *
     * @throws IllegalArgumentException if j is outside [-2, n-1].
     */
    public fun omegaDeriv(j: Int, t: Double): Double {
        requireIndex(j)
        if (t < grid.x(j) || t > grid.x(j + 3)) return 0.0
        val k = intervalOf(t)
        val slot = j - (k - 2)
        if (slot < 0 || slot > 2) return 0.0
        val inv = invM[k]
        val p = frames[k].psiD(t)
        return inv[slot] * p[0] + inv[slot + 3] * p[1] + inv[slot + 6] * p[2]
    }

    /**
     * Second derivative omega_j''(t) (psi_k is replaced by psi_k''). Needed by xi^<0>
     * (de Boor--Fix with r=0). Piecewise constant on the layers; at the grid nodes omega_j'' has a
     * jump (omega_j in C^1 \ C^2), so the value at a node is taken from the right piece.
     *
     * @throws IllegalArgumentException if j is outside [-2, n-1].
     */
    public fun omegaDeriv2(j: Int, t: Double): Double {
        requireIndex(j)
        if (t < grid.x(j) || t > grid.x(j + 3)) return 0.0
        val k = intervalOf(t)
        val slot = j - (k - 2)
        if (slot < 0 || slot > 2) return 0.0
        val inv = invM[k]
        val p = frames[k].psiDD(t)
        return inv[slot] * p[0] + inv[slot + 3] * p[1] + inv[slot + 6] * p[2]
    }

    /**
     * Value of the spline u_h(t) = sum_j c_j omega_j(t), with c of size n+2.
     *
     * @throws IllegalArgumentException if c.size != n + 2 or t is strictly outside the grid interval
     *   (see [intervalOf]).
     */
    public fun evalSpline(c: DoubleArray, t: Double): Double {
        requireCoefficients(c)
        val k = intervalOf(t)
        val w = activeOmega(k, t)
        return c[k] * w[0] + c[k + 1] * w[1] + c[k + 2] * w[2] // indices k-2,k-1,k -> +2
    }

    /**
     * Value of the spline derivative u_h'(t).
     *
     * @throws IllegalArgumentException if c.size != n + 2 or t is strictly outside the grid interval
     *   (see [intervalOf]).
     */
    public fun evalSplineDeriv(c: DoubleArray, t: Double): Double {
        requireCoefficients(c)
        val k = intervalOf(t)
        val inv = invM[k]
        val p = frames[k].psiD(t)
        val w0 = inv[0] * p[0] + inv[3] * p[1] + inv[6] * p[2]
        val w1 = inv[1] * p[0] + inv[4] * p[1] + inv[7] * p[2]
        val w2 = inv[2] * p[0] + inv[5] * p[1] + inv[8] * p[2]
        return c[k] * w0 + c[k + 1] * w1 + c[k + 2] * w2
    }

    /**
     * Value of the second derivative of the spline u_h''(t) (for xi^{<0>} idempotency).
     *
     * @throws IllegalArgumentException if c.size != n + 2 or t is strictly outside the grid interval
     *   (see [intervalOf]).
     */
    public fun evalSplineDeriv2(c: DoubleArray, t: Double): Double {
        requireCoefficients(c)
        val k = intervalOf(t)
        val inv = invM[k]
        val p = frames[k].psiDD(t)
        val w0 = inv[0] * p[0] + inv[3] * p[1] + inv[6] * p[2]
        val w1 = inv[1] * p[0] + inv[4] * p[1] + inv[7] * p[2]
        val w2 = inv[2] * p[0] + inv[5] * p[1] + inv[8] * p[2]
        return c[k] * w0 + c[k + 1] * w1 + c[k + 2] * w2
    }
}

/** The nodes x_j..x_{j+3} are distinct (no multiple nodes merge). */
public fun nonDegenerate(grid: Grid, j: Int): Boolean =
    grid.x(j) < grid.x(j + 1) && grid.x(j + 1) < grid.x(j + 2) && grid.x(j + 2) < grid.x(j + 3)

/** Cross product u x v in R^3 — local arithmetic used to build a_j. */
private fun cross(u: DoubleArray, v: DoubleArray): DoubleArray = doubleArrayOf(
    u[1] * v[2] - u[2] * v[1],
    u[2] * v[0] - u[0] * v[2],
    u[0] * v[1] - u[1] * v[0],
)

/** Dot product in R^3. */
private fun dot(u: DoubleArray, v: DoubleArray): Double = u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

/** All components of the vector are finite (no overflow and no NaN). */
private fun allFinite(u: DoubleArray): Boolean = u.all { it.isFinite() }
