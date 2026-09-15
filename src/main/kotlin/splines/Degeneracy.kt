package splines

import kotlin.math.abs

// ============================================================================
// Degeneracy criterion for scalar denominators
// ============================================================================

/**
 * Relative significance threshold for a scalar quantity (the denominator of the vector `a_j`, the
 * Wronskian, the denominator of `K1`) measured against its own scale.
 *
 * All protected quantities are algebraic sums of products of the input data, so their evaluation is
 * prone to loss of significance through mutual cancellation of the terms. The scale is the sum of the
 * absolute values of the terms (see [cancellationScale]); a quantity is treated as significant when
 * `|value| > DEGENERACY_RELATIVE_EPS * scale` (see [isSignificant]). An absolute threshold is not
 * applicable here: the denominators of the nodal constructions behave like powers of the step `h`,
 * that is, they depend on the scale of the interval.
 *
 * The value 1e-12 separates a quantity lost in round-off noise (of order `1e-15 * scale` for a sum of
 * a few products of two or three factors) from a small but trustworthy quantity. The criterion
 * controls loss of significance, not smallness: a small value computed without cancellation is
 * accepted as significant. As a consequence, on inputs where the terms are large while their
 * difference stays O(1) (for example, divided differences of nodes on intervals with `|x| ≳ 1e12`),
 * the quantity is rejected, since no trustworthy digits are left in it.
 *
 * Approximation-relation matrices `M_k` are not checked by this criterion: their inversion is
 * controlled by a condition number estimate (see [MinimalSplineBasis.MAX_CONDITION]).
 */
internal const val DEGENERACY_RELATIVE_EPS = 1e-12

/**
 * Scale of an algebraic sum: the sum of the absolute values of its terms [terms] — an upper bound for
 * the quantity obtained without mutual cancellation; the round-off noise of the result is
 * proportional to this value.
 */
internal fun cancellationScale(vararg terms: Double): Double {
    var s = 0.0
    for (t in terms) s += abs(t)
    return s
}

/**
 * Significance of the value [value] at the scale [scale]: `|value| > DEGENERACY_RELATIVE_EPS * scale`.
 * The inequality is strict, so for `scale == 0` (all terms are zero) `false` is returned.
 */
internal fun isSignificant(value: Double, scale: Double): Boolean =
    abs(value) > DEGENERACY_RELATIVE_EPS * scale

/** Scale of a dot product in R^3: the sum of the absolute values of the componentwise products. */
internal fun dot3Scale(u: DoubleArray, v: DoubleArray): Double =
    abs(u[0] * v[0]) + abs(u[1] * v[1]) + abs(u[2] * v[2])
