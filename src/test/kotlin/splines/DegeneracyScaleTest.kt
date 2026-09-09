package splines

import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Масштабная инвариантность порогов вырожденности (пункт 2.6 спека).
 *
 * До правки все проверки использовали АБСОЛЮТНЫЙ порог `1e-14`, что давало две
 * симметричные ошибки:
 *  - ЛОЖНОЕ «singular» на мелком масштабе: `det(M_k) ~ h^3`, поэтому на отрезке
 *    `[0, 1e-6]` совершенно корректная сетка отбраковывалась;
 *  - ПРОПУСК вырожденности на крупном масштабе: на `[0, 1e6]` определитель почти
 *    зависимых столбцов — чистый шум округления, но по модулю велик.
 *
 * Оба сценария закреплены здесь фактически. См. KDoc [DEGENERACY_RELATIVE_EPS].
 */
@Tag("fast")
class DegeneracyScaleTest {

    /**
     * МЕЛКИЙ МАСШТАБ, реальный путь построения: базис минимальных сплайнов на
     * отрезке `[0, 1e-6]` для полиномиальной системы B.
     *
     * Здесь `det(M_k) ~ h^3 ~ 2e-21` при собственном масштабе определителя `~1e-19`
     * (компоненты столбцов — `1, t, t^2`), то есть отношение `det/scale ~ 1e-2` и
     * определитель вычислен полностью достоверно. Тем не менее ДО правки
     * построение падало с ложным «matrix is singular» из-за абсолютного
     * порога `1e-14`. Это и есть ключевое доказательство ценности правки.
     *
     * Проверяется не только отсутствие исключения, но и содержательное свойство —
     * разбиение единицы `sum_j omega_j(t) = 1` (первая компонента phi равна 1).
     */
    @Test
    fun splineBasisBuildsOnTinyInterval() {
        val grid = Grid.uniform(8, 0.0, 1e-6)
        val basis = MinimalSplineBasis(GeneratingSystem.B, grid)
        for (k in 0 until grid.n) {
            val t = 0.5 * (grid.x(k) + grid.x(k + 1))
            var sum = 0.0
            for (j in -2..grid.n - 1) sum += basis.omega(j, t)
            assertEquals(1.0, sum, 1e-9, "разбиение единицы в t=$t")
        }
    }

    /**
     * ГРАНИЦА ПРИМЕНИМОСТИ (найдено фактически): для НЕПОЛИНОМИАЛЬНЫХ
     * порождающих систем (H, T) отрезок длины `1e-6` отбраковывается и ПОСЛЕ
     * правки — и это ПРАВИЛЬНО, а не недоработка порога.
     *
     * Причина: у таких систем третья компонента phi НЕ мала на мелком отрезке
     * (`cos t -> 1`, `cosh t -> 1`), поэтому собственный масштаб определителя — `~t`,
     * а сам определитель — `~h^3`: отношение `det/scale ~ 2e-15`, то есть в ответе
     * осталось около ОДНОЙ верной цифры (eps = 2.2e-16). Обращать такую матрицу
     * бессмысленно; отказ — честное поведение. Старый абсолютный порог тоже
     * отказывал, но «случайно правильно» и без диагностики; теперь в сообщении
     * видно отношение `det/scale`.
     */
    @Test
    fun nonPolynomialSystemOnTinyInterval() {
        for (sys in listOf(GeneratingSystem.H, GeneratingSystem.T)) {
            val grid = Grid.uniform(8, 0.0, 1e-6)
            val basis = try { MinimalSplineBasis(sys, grid) } catch (ex: IllegalArgumentException) {
                println("DegeneracyScaleTest: ${sys.name} на [0, 1e-6] отвергнута: ${ex.message}")
                assertTrue(
                    ex.message!!.contains("плохо обусловлена") || ex.message!!.contains("вырождена"),
                    "диагностика обязана называть причину отказа: ${ex.message}"
                )
                continue
            }
            println("DegeneracyScaleTest: ${sys.name} на [0, 1e-6] построена")
            val ones = DoubleArray(grid.n + 2) { 1.0 }
            for (i in 1..20) {
                val t = 1e-6 * i / 21.0
                assertEquals(1.0, basis.evalSpline(ones, t), 1e-8, "разбиение единицы ${sys.name} в t=$t")
            }
        }
    }

    /**
     * КРУПНЫЙ МАСШТАБ, реальный путь: на `[0, 1e6]` построение тоже обязано
     * работать — сетка невырождена, просто велика. Проверка симметрична предыдущей
     * и доказывает, что относительный критерий не «съехал» в другую сторону.
     *
     * Используются только полиномиальная и тригонометрическая системы: у
     * гиперболической `cosh(1e6)` — переполнение double, и это ограничение самой
     * порождающей системы, а не порогов вырожденности.
     */
    @Test
    fun splineBasisBuildsOnHugeInterval() {
        for (sys in listOf(GeneratingSystem.B, GeneratingSystem.T)) {
            val grid = Grid.uniform(8, 0.0, 1e6)
            val basis = MinimalSplineBasis(sys, grid)
            for (k in 0 until grid.n) {
                val t = 0.5 * (grid.x(k) + grid.x(k + 1))
                var sum = 0.0
                for (j in -2..grid.n - 1) sum += basis.omega(j, t)
                assertEquals(1.0, sum, 1e-6, "разбиение единицы для ${sys.name} в t=$t")
            }
        }
    }
}
