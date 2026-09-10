package splines

import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Масштабная инвариантность критериев вырожденности: построение базиса на отрезках длины 1e-6
 * и 1e6. Критерии (относительный порог [DEGENERACY_RELATIVE_EPS] для скалярных знаменателей и
 * число обусловленности матрицы аппроксимационного соотношения в локальных координатах) не
 * зависят от масштаба отрезка: корректная сетка на малом отрезке не отбраковывается, а на большом
 * отрезке шум округления не принимается за значимую величину.
 */
@Tag("fast")
class DegeneracyScaleTest {

    /**
     * Малый масштаб: базис минимальных сплайнов на отрезке `[0, 1e-6]` для полиномиальной
     * системы B строится, и выполнено разбиение единицы `sum_j omega_j(t) = 1` (первая компонента
     * phi равна 1). Абсолютный порог вырожденности отбраковывал бы такую сетку, поскольку
     * `det(M_k) ~ h^3 ~ 2e-21` в глобальных координатах.
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
     * Неполиномиальные системы (H, T) на отрезке длины `1e-6`: построение либо выполняется, и тогда
     * выполнено разбиение единицы, либо отклоняется с диагностикой, называющей причину
     * (вырожденность или плохую обусловленность). Возврат недостоверных значений исключается.
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
     * Большой масштаб: на `[0, 1e6]` построение выполняется, поскольку сетка невырождена.
     *
     * Используются полиномиальная и тригонометрическая системы: для гиперболической `cosh(1e6)`
     * переполняет double, что является ограничением порождающей системы, а не критерия
     * вырожденности.
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
