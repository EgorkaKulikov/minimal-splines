package splines

// ============================================================================
// Grid with multiple boundary nodes
// ============================================================================

/**
 * Grid X: a=x_{-2}=x_{-1}=x_0 < x_1 < ... < x_{n-1} < x_n=x_{n+1}=x_{n+2}=b
 * with nodes of multiplicity 3 at the ends. Stores x_j for j = -2 .. n+2; idx(j)=j+2.
 *
 * @property n number of interior intervals.
 * @param interior interior nodes x_0..x_n (size n+1), including the ends a and b.
 * @throws IllegalArgumentException if n < 1, the size of interior is not n+1, or the nodes are not
 *   strictly increasing.
 */
public class Grid(public val n: Int, interior: DoubleArray) {
    init {
        // For n = 0 there is no interval at all, and the step h = max_j (x_{j+1} - x_j) is undefined.
        require(n >= 1) { "Grid: n must be >= 1 (at least one interval is required), got n=$n" }
        require(interior.size == n + 1) {
            "interior must have size n+1=${n + 1}, got ${interior.size}"
        }
        // Strict increase of the nodes is the invariant behind the step h = max_j (x_{j+1} - x_j),
        // the binary search for an interval in MinimalSplineBasis and the spline supports [x_j, x_{j+3}];
        // it is checked when the grid is built, since the algorithms do not detect its violation.
        for (i in 0 until n) require(interior[i] < interior[i + 1]) {
            "Grid: nodes must be strictly increasing; violated at i=$i: " +
                "x[$i]=${interior[i]} is not less than x[${i + 1}]=${interior[i + 1]}"
        }
    }

    /** Left end of the interval a = x_0. */
    public val a: Double = interior.first()

    /** Right end of the interval b = x_n. */
    public val b: Double = interior.last()

    /**
     * Nodes x_{-2..n+2}, length n+5; idx(j) = j+2. Access only through [x]: the private array
     * preserves the strict-increase invariant checked at construction time.
     */
    private val knots: DoubleArray = DoubleArray(n + 5).also { arr ->
        arr[idx(-2)] = interior[0]; arr[idx(-1)] = interior[0]
        for (j in 0..n) arr[idx(j)] = interior[j]
        arr[idx(n + 1)] = interior[n]; arr[idx(n + 2)] = interior[n]
    }

    /** Step h = max_j (x_{j+1} - x_j) over the interior intervals. */
    public val h: Double = (0 until n).maxOf { interior[it + 1] - interior[it] }

    /**
     * Breakpoint inclusion tolerance: a node is treated as lying strictly inside a subinterval if it is
     * farther than this value from the subinterval's end.
     *
     * Computed as `BREAKPOINT_INCLUSION_EPS_UNIT * max(1, |b - a|)`: the tolerance grows together with
     * the length of the interval, since an absolute tolerance of 1e-15 on an interval of scale 1e3 lies
     * below the ulp of the nodes themselves (`ulp(1e3) ≈ 2.3e-13`) and degenerates into a strict
     * comparison. On intervals shorter than one the tolerance equals `BREAKPOINT_INCLUSION_EPS_UNIT`.
     */
    public val breakpointInclusionEps: Double = BREAKPOINT_INCLUSION_EPS_UNIT * maxOf(1.0, kotlin.math.abs(b - a))

    /** Maps a mathematical index j to an array index. */
    public fun idx(j: Int): Int = j + 2

    /** Node x_j by its mathematical index j (-2..n+2). */
    public fun x(j: Int): Double = knots[idx(j)]

    /**
     * Whether the neighbouring nodes x_j and x_{j+1} coincide, that is, whether the interval
     * [x_j, x_{j+1}] is degenerate.
     *
     * The answer is determined by the index alone, without comparing values: multiple nodes are
     * produced only by the boundary triples x_{-2} = x_{-1} = x_0 = a and x_n = x_{n+1} = x_{n+2} = b,
     * while the interior nodes x_0 < x_1 < ... < x_n increase strictly by the grid invariant. Nodes can
     * coincide only for j < 0 (left multiple node) or j >= n (right multiple node); the comparison
     * `x(j) == x(j+1)` gives the same result.
     *
     * @param j left index of the pair; the admissible range is -2..n+1, for which both nodes x_j and
     *   x_{j+1} lie in the stored range -2..n+2.
     * @throws IllegalArgumentException if j is outside -2..n+1.
     */
    public fun isCoincident(j: Int): Boolean {
        require(j in -2..n + 1) {
            "Grid.isCoincident: j must be in -2..${n + 1} (pair x_j, x_{j+1}), got j=$j"
        }
        return j < 0 || j >= n
    }

    /**
     * Interior nodes x_0..x_n (size n+1) — the smoothness breakpoints of the spline and the boundaries
     * of the quadrature subintervals.
     *
     * The array is not copied; modifying its contents breaks the grid invariants.
     */
    public val breakpoints: DoubleArray = DoubleArray(n + 1) { x(it) }

    /** Grid factories and the base node inclusion tolerance. */
    public companion object {
        /**
         * Base (dimensionless) node inclusion tolerance, referred to the unit interval.
         *
         * The effective tolerance is [breakpointInclusionEps]; this constant is not used in comparisons
         * directly.
         */
        public const val BREAKPOINT_INCLUSION_EPS_UNIT: Double = 1e-15

        /** Uniform grid: x_j = a + (b-a) j/n. */
        public fun uniform(n: Int, a: Double = 0.0, b: Double = 1.0): Grid =
            Grid(n, DoubleArray(n + 1) { a + (b - a) * it / n })

        /**
         * Quasi-uniform grid X^q: x_j = a + (b-a) Psi(j/n), Psi(u) = u + amp sin(2 pi u) — a smooth
         * monotone perturbation of the uniform grid with a local quasi-uniformity parameter that does
         * not depend on n.
         *
         * The function Psi increases strictly only for |amp| < 1/(2 pi) ≈ 0.15915, since
         * Psi'(u) = 1 + 2 pi amp cos(2 pi u). The grid needs a weaker condition — strict increase of the
         * finite set of nodes Psi(j/n): for amp = 1/(2 pi) the nodes increase strictly up to n of about
         * 5·10^5 (beyond that the step near u = 1/2 is lost to rounding), while for amp = 0.16 the nodes
         * increase for n = 8 and stop increasing for n >= 19. For this reason the parameter amp is not
         * restricted separately: strict increase of the constructed nodes is checked by the [Grid]
         * invariant.
         *
         * @param amp perturbation amplitude; any value for which the nodes increase strictly is
         *   admissible.
         * @throws IllegalArgumentException if the constructed nodes are not strictly increasing.
         */
        public fun quasiUniform(n: Int, a: Double = 0.0, b: Double = 1.0, amp: Double = 0.04): Grid =
            Grid(n, DoubleArray(n + 1) { i ->
                val u = i.toDouble() / n
                a + (b - a) * (u + amp * Math.sin(2.0 * Math.PI * u))
            })

        /**
         * Geometric (non-uniform) grid:
         *   x_j = a + (b-a) (q^j - 1)/(q^n - 1),  q = R^{1/(n-1)},  j = 0..n.
         *
         * The ratio of the extreme steps is fixed, h_{n-1}/h_0 = R; the steps grow geometrically
         * (h_j = (b-a) q^j (q-1)/(q^n-1)), so the local quasi-uniformity parameter is constant:
         * h_j/h_{j-1} = q =: mu_n, mu_n -> 1 as n->inf.
         * The grid is substantially non-uniform, yet locally quasi-uniform.
         *
         * @param n number of interior intervals (>= 2).
         * @param R ratio of the extreme steps (> 0). R=1 degenerates into the uniform grid.
         * @throws IllegalArgumentException if n < 2 or R <= 0.
         */
        public fun geometric(n: Int, a: Double = 0.0, b: Double = 1.0, R: Double = 2.0): Grid {
            require(n >= 2) { "geometric: n must be >= 2, got n=$n" }
            require(R > 0.0) { "geometric: R must be > 0, got R=$R" }
            val q = Math.pow(R, 1.0 / (n - 1))
            if (kotlin.math.abs(q - 1.0) < 1e-15) {
                // Degenerate case R=1: the uniform grid (q^n - 1 -> 0).
                return Grid(n, DoubleArray(n + 1) { a + (b - a) * it / n })
            }
            val denom = Math.pow(q, n.toDouble()) - 1.0
            return Grid(n, DoubleArray(n + 1) { j ->
                when (j) {
                    0 -> a
                    n -> b
                    else -> a + (b - a) * (Math.pow(q, j.toDouble()) - 1.0) / denom
                }
            })
        }

        /**
         * Graded two-scale grid with a local quasi-uniformity parameter that does not depend on n. The
         * steps alternate between two values s and ratio*s following the pattern
         * s*(1, ratio, 1, ratio, ...), so the ratio of neighbouring steps h_j/h_{j-1} equals
         * ratio or 1/ratio for every n.
         * Normalisation: the sum of the steps equals b-a (the interval is covered exactly). For odd n
         * the last step remains unpaired (a single value of the pattern).
         *
         * Unlike geometric (mu_n = q -> 1 as n->inf), here mu = ratio is fixed: the grid stays
         * substantially locally non-uniform under refinement (h->0), which is exactly what is needed for
         * a stably large Lebesgue constant. The grid is locally quasi-uniform (neighbouring steps are
         * bounded), so minimal splines on it exist.
         *
         * @param n number of interior intervals (>= 2).
         * @param ratio ratio of neighbouring steps (> 0). ratio=1 degenerates into the uniform grid.
         * @throws IllegalArgumentException if n < 2 or ratio <= 0.
         */
        public fun graded(n: Int, a: Double = 0.0, b: Double = 1.0, ratio: Double = 2.0): Grid {
            require(n >= 2) { "graded: n must be >= 2, got n=$n" }
            require(ratio > 0.0) { "graded: ratio must be > 0, got ratio=$ratio" }
            // Pattern multipliers: even index -> 1, odd index -> ratio (pairs s, ratio*s).
            val mult = DoubleArray(n) { if (it % 2 == 0) 1.0 else ratio }
            val sum = mult.sum()
            val s = (b - a) / sum // base step: normalised to the length of the interval
            val interior = DoubleArray(n + 1)
            interior[0] = a
            var acc = a
            for (i in 0 until n) { acc += s * mult[i]; interior[i + 1] = acc }
            interior[n] = b // guards against accumulated round-off error at the right end
            return Grid(n, interior)
        }
    }
}
