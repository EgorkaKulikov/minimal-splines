package splines

import kotlin.math.cos
import kotlin.math.cosh
import kotlin.math.min
import kotlin.math.sin
import kotlin.math.sinh

// ============================================================================
// 3. ПОРОЖДАЮЩАЯ ВЕКТОР-ФУНКЦИЯ phi(t) = (1, rho(t), sigma(t))^T
// ============================================================================

/**
 * Локальное представление порождающей системы на сеточном интервале: psi(t) = T phi(t), где T —
 * невырожденная матрица 3×3, подобранная так, что компоненты psi и её производных на интервале
 * имеют порядок единицы. Аппроксимационное соотношение sum_j a_j omega_j(t) = phi(t) после левого
 * умножения на T переходит в эквивалентное (T M_k) omega(t) = psi(t) с тем же решением omega(t),
 * поэтому базисные функции от выбора представления не зависят, а обусловленность матрицы T M_k
 * перестаёт зависеть от числа интервалов и от положения и масштаба отрезка.
 *
 * Значения [psi], [psiD], [psiDD] вычисляются непосредственно в локальной переменной, а не как
 * произведение T на phi(t): последнее вернуло бы взаимное сокращение больших слагаемых.
 *
 * @property psi локальная вектор-функция psi(t) = T phi(t);
 * @property psiD её производная psi'(t) = T phi'(t);
 * @property psiDD вторая производная psi''(t) = T phi''(t);
 * @property det определитель матрицы T.
 */
class LocalFrame(
    val psi: (Double) -> DoubleArray,
    val psiD: (Double) -> DoubleArray,
    val psiDD: (Double) -> DoubleArray,
    val det: Double,
)

/**
 * Порождающая вектор-функция phi(t) = (1, rho(t), sigma(t))^T и её производные
 * до второго порядка. phi_0 == 1 обеспечивает разбиение единицы.
 *
 * Необязательный параметр [localFrameFactory] задаёт локальное представление системы на интервале
 * (см. [localFrame]); для встроенных систем [B], [H], [T] он определён.
 */
class GeneratingSystem(
    val name: String,
    val rho: (Double) -> Double,
    val sigma: (Double) -> Double,
    val rhoD: (Double) -> Double,
    val sigmaD: (Double) -> Double,
    val rhoDD: (Double) -> Double,
    val sigmaDD: (Double) -> Double,
    private val localFrameFactory: ((c: Double, h: Double) -> LocalFrame)? = null,
) {
    /** phi(t) = (1, rho(t), sigma(t)). */
    fun phi(t: Double): DoubleArray = doubleArrayOf(1.0, rho(t), sigma(t))

    /** phi'(t) = (0, rho'(t), sigma'(t)). */
    fun phiD(t: Double): DoubleArray = doubleArrayOf(0.0, rhoD(t), sigmaD(t))

    /** phi''(t) = (0, rho''(t), sigma''(t)). */
    fun phiDD(t: Double): DoubleArray = doubleArrayOf(0.0, rhoDD(t), sigmaDD(t))

    /** Вронскиан det(phi, phi', phi'') — проверка невырожденности. */
    fun wronskian(t: Double): Double {
        val a = phi(t); val b = phiD(t); val c = phiDD(t)
        return a[0] * (b[1] * c[2] - b[2] * c[1]) +
            a[1] * (b[2] * c[0] - b[0] * c[2]) +
            a[2] * (b[0] * c[1] - b[1] * c[0])
    }

    /**
     * Локальное представление системы на интервале [c, c + h] (см. [LocalFrame]).
     *
     * Для систем, замкнутых относительно сдвига аргумента (решения линейных однородных
     * дифференциальных уравнений с постоянными коэффициентами, в том числе [B], [H], [T]),
     * представление строится сдвигом в точку c и масштабированием на шаг h. Для системы,
     * не замкнутой относительно сдвига аргумента, локальное представление недоступно:
     * используется глобальное представление с T = I, то есть psi = phi.
     */
    fun localFrame(c: Double, h: Double): LocalFrame =
        localFrameFactory?.invoke(c, h) ?: LocalFrame(::phi, ::phiD, ::phiDD, 1.0)

    companion object {
        /** Полиномиальная phi^B(t) = (1, t, t^2)^T. */
        val B = GeneratingSystem(
            name = "B",
            rho = { t -> t }, sigma = { t -> t * t },
            rhoD = { 1.0 }, sigmaD = { t -> 2.0 * t },
            rhoDD = { 0.0 }, sigmaDD = { 2.0 },
            localFrameFactory = ::polynomialFrame,
        )

        /** Гиперболическая phi^H(t) = (1, sinh t, cosh t)^T. */
        val H = GeneratingSystem(
            name = "H",
            rho = { t -> Math.sinh(t) }, sigma = { t -> Math.cosh(t) },
            rhoD = { t -> Math.cosh(t) }, sigmaD = { t -> Math.sinh(t) },
            rhoDD = { t -> Math.sinh(t) }, sigmaDD = { t -> Math.cosh(t) },
            localFrameFactory = ::hyperbolicFrame,
        )

        /** Тригонометрическая phi^T(t) = (1, sin t, cos t)^T. */
        val T = GeneratingSystem(
            name = "T",
            rho = { t -> Math.sin(t) }, sigma = { t -> Math.cos(t) },
            rhoD = { t -> Math.cos(t) }, sigmaD = { t -> -Math.sin(t) },
            rhoDD = { t -> -Math.sin(t) }, sigmaDD = { t -> -Math.cos(t) },
            localFrameFactory = ::trigonometricFrame,
        )
    }
}

/**
 * Локальное представление системы B: psi(t) = (1, s, s^2), s = (t - c)/h.
 * Матрица T имеет строки (1, 0, 0), (-c/h, 1/h, 0), (c^2/h^2, -2c/h^2, 1/h^2); det T = 1/h^3.
 */
private fun polynomialFrame(c: Double, h: Double): LocalFrame {
    val ih = 1.0 / h
    return LocalFrame(
        psi = { t -> val s = (t - c) * ih; doubleArrayOf(1.0, s, s * s) },
        psiD = { t -> doubleArrayOf(0.0, ih, 2.0 * (t - c) * ih * ih) },
        psiDD = { doubleArrayOf(0.0, 0.0, 2.0 * ih * ih) },
        det = ih * ih * ih,
    )
}

/**
 * Масштаб локального представления для систем H и T: шаг h, но не более единицы — характерного
 * масштаба этих систем; при h > 1 компоненты sinh(t - c), 1 - cos(t - c) уже имеют порядок единицы.
 */
private fun boundedScale(h: Double): Double = min(h, 1.0)

/**
 * Локальное представление системы H: psi(t) = (1, sinh(u)/l, (cosh(u) - 1)/l^2), u = t - c,
 * l = min(h, 1). Третья компонента вычисляется как 2 sinh^2(u/2)/l^2 без потери значимости.
 * Матрица T имеет строки (1, 0, 0), (0, cosh c/l, -sinh c/l), (-1/l^2, -sinh c/l^2, cosh c/l^2);
 * det T = 1/l^3.
 */
private fun hyperbolicFrame(c: Double, h: Double): LocalFrame {
    val il = 1.0 / boundedScale(h)
    val il2 = il * il
    return LocalFrame(
        psi = { t -> val u = t - c; val sh = sinh(0.5 * u); doubleArrayOf(1.0, sinh(u) * il, 2.0 * sh * sh * il2) },
        psiD = { t -> val u = t - c; doubleArrayOf(0.0, cosh(u) * il, sinh(u) * il2) },
        psiDD = { t -> val u = t - c; doubleArrayOf(0.0, sinh(u) * il, cosh(u) * il2) },
        det = il * il2,
    )
}

/**
 * Локальное представление системы T: psi(t) = (1, sin(u)/l, (1 - cos(u))/l^2), u = t - c,
 * l = min(h, 1). Третья компонента вычисляется как 2 sin^2(u/2)/l^2 без потери значимости.
 * Матрица T имеет строки (1, 0, 0), (0, cos c/l, -sin c/l), (1/l^2, -sin c/l^2, -cos c/l^2);
 * det T = -1/l^3.
 */
private fun trigonometricFrame(c: Double, h: Double): LocalFrame {
    val il = 1.0 / boundedScale(h)
    val il2 = il * il
    return LocalFrame(
        psi = { t -> val u = t - c; val s = sin(0.5 * u); doubleArrayOf(1.0, sin(u) * il, 2.0 * s * s * il2) },
        psiD = { t -> val u = t - c; doubleArrayOf(0.0, cos(u) * il, sin(u) * il2) },
        psiDD = { t -> val u = t - c; doubleArrayOf(0.0, -sin(u) * il, cos(u) * il2) },
        det = -il * il2,
    )
}
