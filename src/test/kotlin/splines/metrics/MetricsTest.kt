package splines.metrics

import splines.Grid
import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Тесты метрики ошибки errorEh (равномерная норма на мелкой контрольной сетке).
 * Порядки сходимости `orders` и константа `constCh` живут в `numerical-core`
 * (`numerics.ConvergenceRates`) и проверяются там (`ConvergenceRatesTest`).
 */
@Tag("fast")
class MetricsTest {
    /** errorEh: для exact-eval с известной разностью даёт точный максимум модуля. */
    @Test fun errorEhCapturesMaxDifference() {
        val g = Grid.uniform(4, 0.0, 1.0)
        // eval = exact + bump; максимум |bump| = 0.5 в t=1 (линейный сдвиг 0.5*t)
        val e = errorEh({ t -> t }, { t -> t + 0.5 * t }, g)
        assertEquals(0.5, e, 1e-12)
    }

    /** errorEh = 0 при совпадающих функциях. */
    @Test fun errorEhZeroWhenEqual() {
        val g = Grid.uniform(4)
        assertEquals(0.0, errorEh({ t -> t * t }, { t -> t * t }, g), 1e-15)
    }

    /**
     * errorEh явно требует n >= 1: при n = 0 делитель m = 100n обратился бы в ноль.
     * Конструировать Grid(0, ...) теперь нельзя, поэтому проверяем оба звена контракта:
     * сетка с n = 0 недостижима, а сама метрика считается на минимальной допустимой n = 1.
     */
    @Test fun errorEhRequiresPositiveN() {
        assertFailsWith<IllegalArgumentException> { Grid(0, doubleArrayOf(0.0)) }
        val e = errorEh({ t -> t }, { t -> t + 0.25 }, Grid.uniform(1, 0.0, 1.0))
        assertEquals(0.25, e, 1e-12)
    }

}
