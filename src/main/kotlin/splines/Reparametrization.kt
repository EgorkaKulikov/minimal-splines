package splines

import kotlin.math.expm1
import kotlin.math.ln1p
import kotlin.math.pow

/**
 * Strictly increasing map g on [a, b] used to build the generating system (1, g, g²)
 * (see [GeneratingSystem.reparametrized]); g'(a) = +∞ is allowed.
 *
 * For such a system the minimal splines are omega_j = B_j ∘ g, where B_j are the quadratic B-splines
 * on the knots g(x_j).
 *
 * @property name short name of the map, used in the name of the generating system.
 * @property g the map g(t).
 * @property gD its derivative g'(t); may be +∞ at the left end a.
 * @property gDD its second derivative g''(t); may be infinite at the left end a.
 * @property gInverse the inverse map g^{-1}(y).
 * @property difference g(t) − g(c) computed without cancellation (local-coordinates rule of
 * CONTRIBUTING); the default is the plain difference of the values of g.
 */
public class Reparametrization(
    public val name: String,
    public val g: (Double) -> Double,
    public val gD: (Double) -> Double,
    public val gDD: (Double) -> Double,
    public val gInverse: (Double) -> Double,
    public val difference: (t: Double, c: Double) -> Double = { t, c -> g(t) - g(c) },
) {
    /** Built-in maps. */
    public companion object {
        /**
         * Power map g(t) = (t − a)^beta, 0 < beta ≤ 1, on [a, ∞); for beta < 1 the derivative
         * g'(a) is +∞. The difference g(t) − g(c) is computed as
         * (c − a)^beta · expm1(beta · log1p((t − c)/(c − a))) for c > a and as (t − a)^beta for c = a.
         *
         * @param beta exponent, 0 < beta ≤ 1.
         * @param a left end of the domain (the point where g vanishes).
         * @throws IllegalArgumentException if beta is not in (0, 1] or a is not finite.
         */
        public fun power(beta: Double, a: Double = 0.0): Reparametrization {
            require(beta > 0.0 && beta <= 1.0) { "power: beta must be in (0, 1], got beta=$beta" }
            require(a.isFinite()) { "power: a must be finite, got a=$a" }
            val linear = beta == 1.0
            return Reparametrization(
                name = "pow($beta)",
                g = { t -> (t - a).pow(beta) },
                gD = { t -> if (linear) 1.0 else beta * (t - a).pow(beta - 1.0) },
                gDD = { t -> if (linear) 0.0 else beta * (beta - 1.0) * (t - a).pow(beta - 2.0) },
                gInverse = { y -> a + y.pow(1.0 / beta) },
                difference = { t, c ->
                    val ca = c - a
                    if (ca > 0.0) ca.pow(beta) * expm1(beta * ln1p((t - c) / ca)) else (t - a).pow(beta)
                },
            )
        }
    }
}
