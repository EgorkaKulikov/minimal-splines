package splines

import org.junit.jupiter.api.Tag
import splines.functionals.DeBoorFixFunctionals
import splines.functionals.DiscreteDeBoorFixFunctionals
import splines.metrics.errorEh
import kotlin.math.sinh
import kotlin.test.Test
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Ветви отвержения некорректного входа и численной вырожденности, не достигаемые основными
 * тестами: параметры семейств де Бура--Фикса, параметры фабрик сеток и метрики, вырожденность
 * аппроксимационного соотношения для пользовательской системы в глобальных координатах,
 * переполнение глобального представления H и вырожденность xi^<0> в точке перегиба.
 */
@Tag("fast")
class CoverageGapsTest {
    private val basis = MinimalSplineBasis(GeneratingSystem.B, Grid.uniform(8, 0.0, 1.0))

    private fun message(e: Throwable): String = e.message ?: ""

    @Test
    fun deBoorFix_orderOutsideRangeRejected() {
        for (r in listOf(-1, 3)) {
            val e = assertFailsWith<IllegalArgumentException>("DeBoorFix r=$r") { DeBoorFixFunctionals(basis, r = r) }
            assertTrue(message(e).contains("параметр r должен быть равен 0, 1 или 2") && message(e).contains("получено $r"), message(e))
        }
        for (r in listOf(0, 3)) {
            val e = assertFailsWith<IllegalArgumentException>("DiscreteDeBoorFix r=$r") { DiscreteDeBoorFixFunctionals(basis, r = r) }
            assertTrue(message(e).contains("параметр r должен быть равен 1 или 2") && message(e).contains("получено $r"), message(e))
        }
    }

    @Test
    fun gridFactories_invalidRatioRejected() {
        val geometric = assertFailsWith<IllegalArgumentException> { Grid.geometric(8, R = 0.0) }
        assertTrue(message(geometric).contains("R > 0"), message(geometric))
        val negative = assertFailsWith<IllegalArgumentException> { Grid.geometric(8, R = -2.0) }
        assertTrue(message(negative).contains("R=-2.0"), message(negative))
    }

    @Test
    fun errorEh_invalidRefinementRejected() {
        val e = assertFailsWith<IllegalArgumentException> { errorEh({ 0.0 }, { 0.0 }, basis.grid, refinement = 0) }
        assertTrue(message(e).contains("refinement >= 1"), message(e))
    }

    /**
     * Пользовательская система без локального представления вычисляется в глобальных координатах;
     * для (1, t, t^2) на [0, 1e6] число обусловленности T_k M_k превышает MAX_CONDITION, и конструктор
     * отвергает сетку сообщением о вырожденности матрицы аппроксимационного соотношения.
     */
    @Test
    fun customSystem_illConditionedGlobalMatrixRejected() {
        val globalPolynomial = GeneratingSystem(
            "poly-global", rho = { t -> t }, sigma = { t -> t * t },
            rhoD = { 1.0 }, sigmaD = { t -> 2.0 * t },
            rhoDD = { 0.0 }, sigmaDD = { 2.0 },
        )
        val e = assertFailsWith<IllegalArgumentException> { MinimalSplineBasis(globalPolynomial, Grid.uniform(100, 0.0, 1e6)) }
        assertTrue(
            message(e).contains("Матрица аппроксимационного соотношения на интервале") && message(e).contains("вырождена"),
            message(e),
        )
    }

    /**
     * Базис H на [0, 2000] с h = 1 строится в локальных координатах, но глобальное представление
     * sinh(t), cosh(t) при t > 710 переполняется: диагностический computeA(j) отвергает такой j.
     */
    @Test
    fun globalComputeA_overflowRejected() {
        val hyperbolic = MinimalSplineBasis(GeneratingSystem.H, Grid.uniform(2000, 0.0, 2000.0))
        assertTrue(sinh(1501.0).isInfinite())
        val e = assertFailsWith<IllegalArgumentException> { hyperbolic.computeA(1500) }
        assertTrue(message(e).contains("computeA(j=1500)") && message(e).contains("переполняются"), message(e))
        // Вблизи начала отрезка глобальное представление финитно и совпадает по размеру с локальным.
        assertTrue(hyperbolic.computeA(0).size == 3)
    }

    /**
     * Функционалы xi^<0> требуют rho' sigma'' - rho'' sigma' != 0; для системы (1, t, t^3) эта величина
     * обращается в нуль в точке перегиба t = 0, и при узле сетки в нуле семейство не строится.
     */
    @Test
    fun deBoorFixZero_inflectionNodeRejected() {
        val cubic = GeneratingSystem(
            "cubic", rho = { t -> t }, sigma = { t -> t * t * t },
            rhoD = { 1.0 }, sigmaD = { t -> 3.0 * t * t },
            rhoDD = { 0.0 }, sigmaDD = { t -> 6.0 * t },
        )
        val cubicBasis = MinimalSplineBasis(cubic, Grid.uniform(8, -1.0, 1.0))
        val e = assertFailsWith<IllegalArgumentException> { DeBoorFixFunctionals(cubicBasis, r = 0) }
        assertTrue(message(e).contains("вырожденный знаменатель rho'sigma''-rho''sigma'"), message(e))
    }
}
