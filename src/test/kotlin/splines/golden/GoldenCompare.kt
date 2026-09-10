package splines.golden

import splines.golden.GoldenIo.assertBits
import splines.golden.GoldenIo.assertClose
import splines.golden.GoldenIo.toDbl
import splines.golden.GoldenIo.toVec
import kotlin.test.assertEquals
import kotlin.test.assertNull
import kotlin.test.fail

/** Рекурсивное сравнение деревьев эталона и повторного вычисления. */
object GoldenCompare {
    enum class Mode { BITS, REL }

    const val REL_TOL = 1e-12

    /**
     * Относительный допуск для наборов `basis` и `functionals`.
     *
     * Эталоны 0.1.0 получены обращением M_k по формулам Крамера в глобальных координатах
     * порождающей системы при числах обусловленности до ~10⁴ (сетки с n ≤ 32, включая
     * градуированные с мелкими интервалами); погрешность такого пути — порядка cond·ε ≈ 10⁻¹²…10⁻¹¹,
     * что и наблюдается при сравнении с LAPACK (макс. 1.05·10⁻¹¹). Допуск 10⁻¹⁰ покрывает
     * погрешность эталона с десятикратным запасом и при этом на порядки меньше любой
     * погрешности метода.
     */
    const val BASIS_TOLERANCE = 1e-10

    /**
     * Относительный допуск для величин, вычисляемых через `Math.sin`, `Math.cos`, `Math.pow`
     * (порождающая система T, сетки `quasiUniform` и `geometric`): HotSpot использует
     * платформенные интринсики, результаты которых на x86_64 и aarch64 различаются в младшем
     * бите. В сетке `geometric` расхождение в одну единицу младшего разряда у q = R^{1/(n-1)}
     * переносится на q^j с множителем j ≤ n, а расхождения отдельных степеней складываются:
     * при n = 32 это до ~7 единиц младшего разряда единицы. Допуск в 16 единиц младшего
     * разряда покрывает это различие с запасом и исключает любое расхождение алгоритма.
     */
    val TRANSCENDENTAL_TOLERANCE: Double = 16 * Math.ulp(1.0)

    private fun isHex(s: String): Boolean = s.length == 16 && s.all { it in '0'..'9' || it in 'a'..'f' }

    fun compare(expected: Any?, got: Any?, label: String, mode: Mode, tol: Double = REL_TOL) {
        when (expected) {
            null -> assertNull(got, label)
            is Map<*, *> -> {
                val g = got as? Map<*, *> ?: fail("$label: ожидался объект, получено $got")
                assertEquals(expected.keys.map { it.toString() }.toSet(), g.keys.map { it.toString() }.toSet(), "$label: набор ключей")
                for ((k, v) in expected) compare(v, g[k], "$label.$k", mode, tol)
            }
            is List<*> -> {
                val g = got as? List<*> ?: fail("$label: ожидался список, получено $got")
                assertEquals(expected.size, g.size, "$label: размер")
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
