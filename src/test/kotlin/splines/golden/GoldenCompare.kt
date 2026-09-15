package splines.golden

import splines.golden.GoldenIo.assertBits
import splines.golden.GoldenIo.assertClose
import splines.golden.GoldenIo.toDbl
import splines.golden.GoldenIo.toVec
import kotlin.test.assertEquals
import kotlin.test.assertNull
import kotlin.test.fail

/** Recursive comparison of the golden reference tree against a recomputed one. */
object GoldenCompare {
    enum class Mode { BITS, REL }

    const val REL_TOL = 1e-12

    /**
     * Relative tolerance for the `basis` and `functionals` sets.
     *
     * The 0.1.0 golden references were produced by inverting M_k with Cramer's rule in the global
     * coordinates of the generating system, at condition numbers up to ~10⁴ (grids with n ≤ 32,
     * including graded ones with small intervals); the error of that path is of order cond·ε ≈ 10⁻¹²…10⁻¹¹,
     * which is exactly what the comparison with LAPACK shows (max. 1.05·10⁻¹¹). A tolerance of 10⁻¹⁰ covers
     * the golden reference error with a tenfold margin while still being orders of magnitude smaller than
     * any method error.
     */
    const val BASIS_TOLERANCE = 1e-10

    /**
     * Relative tolerance for quantities computed through `Math.sin`, `Math.cos`, `Math.pow`
     * (generating system T, grids `quasiUniform` and `geometric`): HotSpot uses
     * platform intrinsics whose results differ in the last bit between x86_64 and aarch64.
     * In the `geometric` grid a one-ulp discrepancy in q = R^{1/(n-1)} is carried over to q^j with a
     * factor of j ≤ n, and the discrepancies of the individual powers accumulate:
     * at n = 32 this amounts to ~7 ulps of unity. A tolerance of 16 ulps
     * covers this difference with a margin and rules out any algorithmic discrepancy.
     */
    val TRANSCENDENTAL_TOLERANCE: Double = 16 * Math.ulp(1.0)

    private fun isHex(s: String): Boolean = s.length == 16 && s.all { it in '0'..'9' || it in 'a'..'f' }

    fun compare(expected: Any?, got: Any?, label: String, mode: Mode, tol: Double = REL_TOL) {
        when (expected) {
            null -> assertNull(got, label)
            is Map<*, *> -> {
                val g = got as? Map<*, *> ?: fail("$label: expected an object, got $got")
                assertEquals(expected.keys.map { it.toString() }.toSet(), g.keys.map { it.toString() }.toSet(), "$label: key set")
                for ((k, v) in expected) compare(v, g[k], "$label.$k", mode, tol)
            }
            is List<*> -> {
                val g = got as? List<*> ?: fail("$label: expected a list, got $got")
                assertEquals(expected.size, g.size, "$label: size")
                val hexList = expected.isNotEmpty() && expected.all { it is String && isHex(it) } &&
                    g.all { it is String && isHex(it) }
                if (hexList) {
                    val e = toVec(expected); val gg = toVec(g)
                    if (mode == Mode.BITS) for (i in e.indices) assertBits(e[i], gg[i], "$label[$i]")
                    else assertClose(e, gg, tol, label)
                } else {
                    for (i in expected.indices) compare(expected[i], g[i], "$label[$i]", mode, tol)
                }
            }
            is String -> {
                if (isHex(expected) && got is String && isHex(got)) {
                    if (mode == Mode.BITS) assertBits(toDbl(expected), toDbl(got), label)
                    else assertClose(toDbl(expected), toDbl(got), tol, label)
                } else assertEquals(expected, got, label)
            }
            else -> assertEquals(expected, got, label)
        }
    }
}
