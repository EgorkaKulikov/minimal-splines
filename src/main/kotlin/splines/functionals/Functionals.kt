package splines.functionals

import kotlin.math.abs

import numerics.DenseMatrix
import numerics.LinearAlgebra
import numerics.NumericsContext
import splines.DEGENERACY_RELATIVE_EPS
import splines.Grid
import splines.MinimalSplineBasis
import splines.cancellationScale
import splines.isSignificant

// ============================================================================
// Семейства (квази)проекционных функционалов: theta (проекционные), xi (де Бура–Фикса,
// значение и производные), xitilde (дискретизованные де Бура–Фикса), mu (усредняющие),
// lambda (трёхточечные). Источники формул указаны в заголовках разделов.
// ============================================================================

/**
 * Аппроксимационный функционал chi_j. Семейства theta, mu, lambda используют только значения f;
 * семейство xi использует также производные f' и f'', поэтому интерфейс принимает функцию вместе
 * с производными, а семейства без производных их не используют.
 */
public interface ApproxFunctional {
    /**
     * Значение chi_j(f) по функции f, её первой производной fD и второй производной fDD.
     * Семейства theta, mu, lambda используют только f; xi^<1>, xi^<2> — также fD; xi^<0> — также fDD.
     * По умолчанию fDD равна нулю, что достаточно для семейств, не использующих вторую производную.
     */
    public fun apply(f: (Double) -> Double, fD: (Double) -> Double, fDD: (Double) -> Double = { 0.0 }): Double

    /**
     * Сумма модулей коэффициентов функционала в его представлении, без нормировки на шаг сетки.
     *
     * Для функционалов-значений (`chi_j(f) = sum_k c_k f(t_k)`; семейства theta, mu, lambda,
     * xitilde) величина совпадает с нормой функционала относительно возмущений входных значений:
     * возмущение данных величины eps в норме максимума даёт погрешность не более `absSum() * eps`.
     * Для функционалов с производными ([DerivFunctional], [SecondDerivFunctional]) величина
     * нормой относительно возмущений данных не является (см. KDoc этих классов) и служит только
     * диагностикой представления.
     */
    public fun absSum(): Double
}

/** Удобная обёртка: chi_j(f) без явных производных (производные = 0). */
public fun ApproxFunctional.apply(f: (Double) -> Double): Double = apply(f, { 0.0 }, { 0.0 })

/**
 * Семейство аппроксимационных функционалов {chi_j}_{j=-2}^{n-1} и (квази)проектор
 * P_chi g = sum_j chi_j(g) omega_j. Общий интерфейс семейств theta, xi, xitilde, mu, lambda.
 *
 * @property basis базис минимальных сплайнов, к которому относятся функционалы.
 * @property name имя семейства (theta, xi, xi<r>, xitilde, xitilde<r>, mu, lambda).
 * @property isProjector `true` для проекторов (theta, xi): выполнена биортогональность
 *   chi_i(omega_j) = delta_ij и P_chi^2 = P_chi; `false` для квазиинтерполянтов (xitilde, mu, lambda).
 * @property usesDerivative `true` для семейства xi, использующего производную образа.
 */
public abstract class FunctionalFamily(
    public val basis: MinimalSplineBasis,
    public val name: String,
    /**
     * Контекст численных вычислений: реализация BLAS/LAPACK, которой семейства theta, mu, lambda
     * решают в конструкторе СЛАУ 3×3 и 5×5. Параметр принимают все семейства, в том числе не
     * использующие линейную алгебру (xi, xitilde), что допускает единообразную сверку контекста
     * любого семейства с контекстом вызывающего кода.
     */
    public val ctx: NumericsContext = NumericsContext.default(),
) {
    /** Сетка базиса. */
    public val grid: Grid = basis.grid

    /** Число интервалов сетки; семейство состоит из n + 2 функционалов. */
    public val n: Int = grid.n
    public abstract val isProjector: Boolean
    public abstract val usesDerivative: Boolean

    /** `true` для семейств, использующих вторую производную образа (xi^<0>). */
    public open val usesSecondDerivative: Boolean = false

    /** Функционал chi_j, j = -2..n-1. */
    public abstract fun chi(j: Int): ApproxFunctional

    /** Коэффициенты проекции P_chi g = sum chi_j(g) omega_j: вектор (chi_j(g)) размера n+2. */
    public fun projectorCoeffs(
        g: (Double) -> Double,
        gD: (Double) -> Double = { 0.0 },
        gDD: (Double) -> Double = { 0.0 },
    ): DoubleArray = DoubleArray(n + 2) { chi(it - 2).apply(g, gD, gDD) }

    /**
     * Максимум [ApproxFunctional.absSum] по всем j = -2..n-1.
     *
     * Для семейств из функционалов-значений (`usesDerivative == false`) это константа `C_chi` —
     * оценка усиления возмущения входных данных (квази)проектором P_chi в норме максимума.
     * Для семейства xi (`usesDerivative == true`) величина такой оценкой не является: коэффициенты
     * при производных убывают с шагом сетки, тогда как усиление возмущения при численном
     * дифференцировании растёт (см. KDoc [DerivFunctional]).
     */
    public fun cChi(): Double = (-2..n - 1).maxOf { chi(it).absSum() }
}

// ----------------------------------------------------------------------------
// theta — проекционные функционалы
// Источник: Kulikov, Makarov (Записки научных семинаров ПОМИ, 2025, т. 542,
// с. 126–143). См. docs/REFERENCES.md, раздел 2.
// ----------------------------------------------------------------------------

/**
 * Функционал-значение: линейная комбинация значений f в точках [nodes] с коэффициентами [coeffs].
 *
 * Массивы не копируются; изменение их содержимого нарушает инварианты функционала.
 *
 * @property nodes опорные точки.
 * @property coeffs коэффициенты при значениях в [nodes].
 * @throws IllegalArgumentException если длины массивов различны.
 */
public class ValueFunctional(public val nodes: DoubleArray, public val coeffs: DoubleArray) : ApproxFunctional {
    init { require(nodes.size == coeffs.size) }
    override fun apply(f: (Double) -> Double, fD: (Double) -> Double, fDD: (Double) -> Double): Double {
        var s = 0.0
        for (k in nodes.indices) s += coeffs[k] * f(nodes[k])
        return s
    }
    /**
     * Сумма модулей [coeffs]. Для функционала-значения это в точности норма
     * функционала как коэффициента усиления возмущения входных значений:
     * `|chi(f + e) - chi(f)| <= absSum() * max|e|`.
     */
    override fun absSum(): Double = coeffs.fold(0.0) { acc, v -> acc + abs(v) }
}

/**
 * Семейство проекционных функционалов theta_j.
 *
 * Внутренние и краевые функционалы строятся локальной биортогонализацией: решается
 * система theta_j(omega_i) = delta_ij по узлам сетки и серединам интервалов. Это
 * устойчивое эквивалентное представление закрытой формулы из источника; совпадение
 * двух представлений проверяется тестом (см. [closedFormInternal]).
 *
 * Семейство является проектором: выполнена биортогональность theta_i(omega_j) = delta_ij,
 * откуда P_theta^2 = P_theta. Краевые функционалы (j = -2 и j = n-1) — значения f(x_0) и f(x_n)
 * согласно определению в источнике.
 *
 * @param ctx контекст численных вычислений для СЛАУ локальной биортогонализации.
 */
public class ProjFunctionals(
    basis: MinimalSplineBasis,
    ctx: NumericsContext = NumericsContext.default(),
) : FunctionalFamily(basis, "theta", ctx) {
    override val isProjector: Boolean = true
    override val usesDerivative: Boolean = false
    private val funcs: Array<ApproxFunctional> = Array(n + 2) { buildTheta(it - 2) }
    override fun chi(j: Int): ApproxFunctional = funcs[j + 2]

    private fun mid(p: Double, q: Double) = 0.5 * (p + q)

    private fun buildTheta(j: Int): ApproxFunctional {
        val x0 = grid.x(0); val x1 = grid.x(1)
        val xnm1 = grid.x(n - 1); val xn = grid.x(n)
        return when (j) {
            -2 -> ValueFunctional(doubleArrayOf(x0), doubleArrayOf(1.0))
            -1 -> localFunctional(j = -1, points = doubleArrayOf(x0, mid(x0, x1), x1), indices = intArrayOf(-2, -1, 0))
            n - 1 -> ValueFunctional(doubleArrayOf(xn), doubleArrayOf(1.0))
            n - 2 -> localFunctional(j = n - 2, points = doubleArrayOf(xnm1, mid(xnm1, xn), xn), indices = intArrayOf(n - 3, n - 2, n - 1))
            else -> {
                val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2); val xj3 = grid.x(j + 3)
                localFunctional(
                    j = j,
                    points = doubleArrayOf(xj, mid(xj, xj1), mid(xj1, xj2), mid(xj2, xj3), xj3),
                    indices = intArrayOf(j - 2, j - 1, j, j + 1, j + 2),
                )
            }
        }
    }

    /** Локальная биортогонализация: coeff так, что sum_p coeff_p omega_i(points_p)=delta_ij. */
    private fun localFunctional(j: Int, points: DoubleArray, indices: IntArray): ValueFunctional {
        val m = points.size
        val matrix = DenseMatrix.build(m, m) { r, c -> basis.omega(indices[r], points[c]) }
        val rhs = DoubleArray(m) { if (indices[it] == j) 1.0 else 0.0 }
        val coeff = LinearAlgebra.solve(matrix, rhs, ctx.backend)
        return ValueFunctional(points, coeff)
    }

    /**
     * Закрытая (явная) формула функционала theta_j для внутреннего индекса j.
     *
     * Используется только как независимая сверка с основным построением через
     * локальную биортогонализацию: оба представления обязаны совпадать.
     *
     * Функционал опирается на пять точек: узлы x_j, x_{j+3} и три середины
     * интервалов носителя. Обозначения соответствуют источнику (docs/REFERENCES.md,
     * раздел 2): значения omega_j в этих точках образуют величины A..E, а знаменатель
     * K1 = C^2 D - B C E - A D E - D E^2.
     *
     * @param j внутренний индекс функционала (краевые j не поддерживаются: у них иная ветвь формулы)
     * @throws IllegalStateException если знаменатель вырождается
     */
    internal fun closedFormInternal(j: Int): ValueFunctional {
        val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2); val xj3 = grid.x(j + 3)
        val a = basis.omega(j, mid(xj, xj1))
        val c = basis.omega(j, mid(xj1, xj2))
        val e = basis.omega(j, mid(xj2, xj3))
        val b = basis.omega(j, xj1)
        val d = basis.omega(j, xj2)
        val k1 = c * c * d - b * c * e - a * d * e - d * e * e
        // Масштабом знаменателя служит сумма модулей четырёх слагаемых K1: значения omega_j имеют
        // порядок единицы, а их разности малы как степени h, поэтому абсолютный порог не инвариантен
        // к шагу сетки. См. KDoc splines.DEGENERACY_RELATIVE_EPS.
        val k1Scale = cancellationScale(c * c * d, b * c * e, a * d * e, d * e * e)
        check(isSignificant(k1, k1Scale)) {
            "closedFormInternal(j=$j): вырожденный знаменатель K1=$k1, scale=$k1Scale " +
                "(значимость потеряна: порог $DEGENERACY_RELATIVE_EPS)"
        }
        return ValueFunctional(
            doubleArrayOf(xj, mid(xj, xj1), mid(xj1, xj2), mid(xj2, xj3), xj3),
            doubleArrayOf(e * e / k1, -d * e / k1, (c * d - b * e) / k1, -d * e / k1, e * e / k1),
        )
    }
}
// ----------------------------------------------------------------------------
// xi — функционалы де Бура–Фикса (значение и производные)
// Источник: Kulikov, Makarov, On de Boor–Fix Type Functionals for Minimal Splines
// (Topics in Classical and Modern Analysis, Springer, 2019, p. 211–225).
// См. docs/REFERENCES.md, раздел 2.
// ----------------------------------------------------------------------------

/**
 * Функционал вида xi(u) = u(node) + cD u'(node) (де Бура–Фикса, r = 1, 2).
 * При cD = 0 сводится к значению u(node) — краевой функционал u(x_0), u(x_n).
 *
 * @property node узел функционала.
 * @property cD коэффициент при производной; имеет размерность длины и порядок шага сетки.
 */
public class DerivFunctional(public val node: Double, public val cD: Double) : ApproxFunctional {
    override fun apply(f: (Double) -> Double, fD: (Double) -> Double, fDD: (Double) -> Double): Double =
        f(node) + cD * fD(node)

    /**
     * Сумма модулей коэффициентов представления `1 + |cD|`.
     *
     * Величина не оценивает усиление возмущения данных: коэффициент `cD` имеет порядок шага `h`,
     * поэтому `1 + |cD| -> 1` при `h -> 0`, тогда как возмущение eps, пропущенное через численное
     * дифференцирование, усиливается как `eps/h`. Величина пригодна для диагностики представления
     * (порядок коэффициента при производной), но не для оценок влияния шума данных.
     */
    override fun absSum(): Double = 1.0 + abs(cD)
}

/**
 * Функционал вида xi^<0>(u) = u(node) + c1 u'(node) + c2 u''(node) (де Бура–Фикса, r = 0):
 * использует значение, первую и вторую производные в одном узле.
 *
 * @property node узел функционала.
 * @property c1 коэффициент при первой производной (порядок h).
 * @property c2 коэффициент при второй производной (порядок h^2).
 */
public class SecondDerivFunctional(public val node: Double, public val c1: Double, public val c2: Double) : ApproxFunctional {
    override fun apply(f: (Double) -> Double, fD: (Double) -> Double, fDD: (Double) -> Double): Double =
        f(node) + c1 * fD(node) + c2 * fDD(node)

    /**
     * Сумма модулей коэффициентов представления `1 + |c1| + |c2|`.
     *
     * Величина не оценивает усиление возмущения данных по той же причине, что и в
     * [DerivFunctional]: `c1` имеет порядок `h`, `c2` — порядок `h^2`, и сумма стремится к 1 при
     * измельчении сетки, тогда как усиление возмущения первой и второй производными растёт как
     * `1/h` и `1/h^2`. Пригодна только для диагностики представления.
     */
    override fun absSum(): Double = 1.0 + abs(c1) + abs(c2)
}

/**
 * Семейство функционалов де Бура–Фикса xi_j^{<r>}, r из {0, 1, 2}.
 *
 * Все три семейства — проекторы: выполнена биортогональность xi_i(omega_j) = delta_ij.
 *
 *  - xi^<1>(u) = u(x_{j+1}) + C1_j u'(x_{j+1}),
 *      C1_j = ((sigma_{j+2}-sigma_{j+1})rho'_{j+2} - (rho_{j+2}-rho_{j+1})sigma'_{j+2}) / W_j;
 *  - xi^<2>(u) = u(x_{j+2}) + C2_j u'(x_{j+2}),
 *      C2_j = ((sigma_{j+2}-sigma_{j+1})rho'_{j+1} - (rho_{j+2}-rho_{j+1})sigma'_{j+1}) / W_j;
 *      где W_j = rho'_{j+2}sigma'_{j+1} - rho'_{j+1}sigma'_{j+2} — вронскиан;
 *  - xi^<0>(u) = u(x_j) + (N1_j/Delta_j) u'(x_j) + (N2_j/Delta_j) u''(x_j)
 *      — использует вторую производную образа, то есть требует u из C^2.
 *
 * Коэффициенты вычисляются в локальных координатах интервала (x_{j+1}, x_{j+2}) (см. [splines.LocalFrame]):
 * формулы для rho, sigma применимы к компонентам psi_1, psi_2, поскольку первая строка матрицы T
 * равна (1, 0, 0) и psi порождает то же пространство. Это исключает потерю значимости при малой
 * длине отрезка или его удалении от нуля; результат от выбора координат не зависит, так как
 * функционал определён условием биортогональности к базису.
 *
 * Краевые функционалы j = -2 и j = n-1 — значения u(x_0) и u(x_n); биортогональность
 * при этом выборе выполнена для всех r и всех порождающих систем.
 *
 * @property r порядок функционала из {0, 1, 2}; по умолчанию r = 1.
 * @param ctx контекст численных вычислений (линейная алгебра семейством не используется).
 * @throws IllegalArgumentException если r вне {0, 1, 2} либо на некотором интервале вронскиан
 *   или знаменатель Delta_j вырожден (см. [splines.DEGENERACY_RELATIVE_EPS]).
 */
public class DeBoorFixFunctionals(
    basis: MinimalSplineBasis,
    public val r: Int = 1,
    ctx: NumericsContext = NumericsContext.default(),
) : FunctionalFamily(basis, if (r == 1) "xi" else "xi<$r>", ctx) {
    init { require(r in 0..2) { "DeBoorFix: r must be in {0,1,2}, got $r" } }
    override val isProjector: Boolean = true
    override val usesDerivative: Boolean = true
    override val usesSecondDerivative: Boolean = (r == 0)
    private val funcs: Array<ApproxFunctional> = Array(n + 2) { buildXi(it - 2) }
    override fun chi(j: Int): ApproxFunctional = funcs[j + 2]

    private fun buildXi(j: Int): ApproxFunctional {
        if (j == -2) return DerivFunctional(grid.x(0), 0.0)
        if (j == n - 1) return DerivFunctional(grid.x(n), 0.0)
        return when (r) {
            0 -> buildXi0(j)
            2 -> buildXi2(j)
            else -> buildXi1(j)
        }
    }

    /** xi^<1>_j: узел x_{j+1}, коэффициент C1_j. */
    private fun buildXi1(j: Int): ApproxFunctional {
        val x1 = grid.x(j + 1); val x2 = grid.x(j + 2)
        // Компоненты psi_1, psi_2 локальной системы интервала (x_{j+1}, x_{j+2}) и их
        // производные: формула та же, что для rho, sigma, поскольку psi = T phi с первой
        // строкой T, равной (1, 0, 0), порождает то же пространство.
        val fr = basis.frame(j + 1)
        val p1 = fr.psi(x1); val p2 = fr.psi(x2)
        val d1 = fr.psiD(x1); val d2 = fr.psiD(x2)
        val rho1 = p1[1]; val rho2 = p2[1]
        val sig1 = p1[2]; val sig2 = p2[2]
        val rhoD1 = d1[1]; val rhoD2 = d2[1]
        val sigD1 = d1[2]; val sigD2 = d2[2]
        // Вронскиан W_j — разность двух произведений; масштаб = сумма их модулей.
        // Так проверка не зависит от масштаба самих rho', sigma' (для системы H они
        // растут как cosh, для T ограничены единицей, а на мелкой сетке разность мала).
        val denom = rhoD2 * sigD1 - rhoD1 * sigD2
        val denomScale = cancellationScale(rhoD2 * sigD1, rhoD1 * sigD2)
        require(isSignificant(denom, denomScale)) {
            "buildXi1(j=$j): degenerate Wronskian rhoD2*sigD1 - rhoD1*sigD2=$denom, " +
                "scale=$denomScale (значимость потеряна: порог $DEGENERACY_RELATIVE_EPS)"
        }
        val cD = ((sig2 - sig1) * rhoD2 - (rho2 - rho1) * sigD2) / denom
        return DerivFunctional(x1, cD)
    }

    /** xi^<2>_j: узел x_{j+2}, коэффициент C2_j (тот же знаменатель W_j, штрихи в x_{j+1}). */
    private fun buildXi2(j: Int): ApproxFunctional {
        val x1 = grid.x(j + 1); val x2 = grid.x(j + 2)
        // Компоненты psi_1, psi_2 локальной системы интервала (x_{j+1}, x_{j+2}) и их
        // производные: формула та же, что для rho, sigma, поскольку psi = T phi с первой
        // строкой T, равной (1, 0, 0), порождает то же пространство.
        val fr = basis.frame(j + 1)
        val p1 = fr.psi(x1); val p2 = fr.psi(x2)
        val d1 = fr.psiD(x1); val d2 = fr.psiD(x2)
        val rho1 = p1[1]; val rho2 = p2[1]
        val sig1 = p1[2]; val sig2 = p2[2]
        val rhoD1 = d1[1]; val rhoD2 = d2[1]
        val sigD1 = d1[2]; val sigD2 = d2[2]
        // Тот же вронскиан, что и в buildXi1: масштаб — сумма модулей двух произведений.
        val denom = rhoD2 * sigD1 - rhoD1 * sigD2
        val denomScale = cancellationScale(rhoD2 * sigD1, rhoD1 * sigD2)
        require(isSignificant(denom, denomScale)) {
            "buildXi2(j=$j): degenerate Wronskian rhoD2*sigD1 - rhoD1*sigD2=$denom, " +
                "scale=$denomScale (значимость потеряна: порог $DEGENERACY_RELATIVE_EPS)"
        }
        val cD = ((sig2 - sig1) * rhoD1 - (rho2 - rho1) * sigD1) / denom
        return DerivFunctional(x2, cD)
    }

    /**
     * xi^<0>_j: узел x_j, коэффициенты N1_j/Delta_j (при u') и N2_j/Delta_j (при u'') по формулам
     * источника; Delta_j = W_j (rho'_j sigma''_j - rho''_j sigma'_j).
     */
    private fun buildXi0(j: Int): ApproxFunctional {
        val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2)
        // Значения psi_1, psi_2 и производные в трёх узлах в локальных координатах
        // интервала (x_{j+1}, x_{j+2}); узел x_j отстоит от начала отсчёта на h_j.
        val fr = basis.frame(j + 1)
        val pj = fr.psi(xj); val dj = fr.psiD(xj); val ddj = fr.psiDD(xj)
        val pj1 = fr.psi(xj1); val dj1 = fr.psiD(xj1)
        val pj2 = fr.psi(xj2); val dj2 = fr.psiD(xj2)
        val rj = pj[1]; val sj = pj[2]
        val rDj = dj[1]; val sDj = dj[2]
        val rDDj = ddj[1]; val sDDj = ddj[2]
        val rj1 = pj1[1]; val sj1 = pj1[2]
        val rDj1 = dj1[1]; val sDj1 = dj1[2]
        val rj2 = pj2[1]; val sj2 = pj2[2]
        val rDj2 = dj2[1]; val sDj2 = dj2[2]

        // Delta_j — произведение двух миноров 2x2 и вырождается ровно тогда, когда вырожден один
        // из множителей. Множители проверяются по отдельности, каждый на своём масштабе: так
        // диагностика указывает причину, а критерий не зависит от произведения масштабов.
        val wronskian12 = rDj1 * sDj2 - rDj2 * sDj1
        val wronskian12Scale = cancellationScale(rDj1 * sDj2, rDj2 * sDj1)
        require(isSignificant(wronskian12, wronskian12Scale)) {
            "buildXi0(j=$j): degenerate Wronskian W_j=$wronskian12, scale=$wronskian12Scale " +
                "(значимость потеряна: порог $DEGENERACY_RELATIVE_EPS)"
        }
        val curvature = rDj * sDDj - rDDj * sDj
        val curvatureScale = cancellationScale(rDj * sDDj, rDDj * sDj)
        require(isSignificant(curvature, curvatureScale)) {
            "buildXi0(j=$j): degenerate rho'sigma''-rho''sigma'=$curvature, scale=$curvatureScale " +
                "(значимость потеряна: порог $DEGENERACY_RELATIVE_EPS)"
        }
        val delta = wronskian12 * curvature
        // Проверки множителей не контролируют само произведение: два значимых, но малых множителя
        // (порядка 1e-200) дают delta == 0.0 при исчезновении порядка, два больших — Inf при
        // переполнении; в обоих случаях n1/delta оказалось бы нечисловым. Используется `require`,
        // как и для множителей: все случаи негодного входа завершаются IllegalArgumentException.
        require(delta.isFinite() && delta != 0.0) {
            "buildXi0(j=$j): Delta_j=$delta непригодна как знаменатель (underflow/overflow произведения " +
                "значимых по отдельности множителей): W_j=$wronskian12, " +
                "rho'sigma''-rho''sigma'=$curvature"
        }
        val n1 = (rj1 * sDj1 - rDj1 * sj1) * (rDDj * sDj2 - rDj2 * sDDj) +
            (rDj1 * sDj2 - rDj2 * sDj1) * (rDDj * sj - rj * sDDj) +
            (rj2 * sDj2 - rDj2 * sj2) * (rDj1 * sDDj - rDDj * sDj1)
        val n2 = (rj1 * sDj1 - rDj1 * sj1) * (rDj2 * sDj - rDj * sDj2) +
            (rDj1 * sDj2 - rDj2 * sDj1) * (rj * sDj - rDj * sj) +
            (rj2 * sDj2 - rDj2 * sj2) * (rDj * sDj1 - rDj1 * sDj)
        // Конечности знаменателя недостаточно: при субнормальном ненулевом delta (порядка 1e-310)
        // частные n1/delta, n2/delta переполняются, поэтому проверяется сам результат.
        val c1 = n1 / delta
        val c2 = n2 / delta
        require(c1.isFinite() && c2.isFinite()) {
            "buildXi0(j=$j): коэффициенты нечисловые: N1/Delta=$c1, N2/Delta=$c2 " +
                "(N1=$n1, N2=$n2, Delta_j=$delta — переполнение при делении на субнормальный знаменатель)"
        }
        return SecondDerivFunctional(xj, c1, c2)
    }
}

// ----------------------------------------------------------------------------
// xitilde — дискретизованные функционалы де Бура–Фикса (без производных)
// ----------------------------------------------------------------------------

/**
 * Дискретизованные функционалы де Бура–Фикса xitilde^{<r>}_j, r из {1, 2}, не использующие
 * производную.
 *
 * Производная f'(x_k) в xi^{<1>}, xi^{<2>} заменена центральной разделённой разностью узловых
 * значений f'(x_k) ≈ (f(x_{k+1}) - f(x_{k-1})) / (x_{k+1} - x_{k-1}), корректной и на
 * неравномерных сетках:
 *   xitilde^{<1>}_j(f) = f(x_{j+1}) + w1_j (f(x_{j+2}) - f(x_j)) / (x_{j+2} - x_j),
 *   xitilde^{<2>}_j(f) = f(x_{j+2}) + w2_j (f(x_{j+3}) - f(x_{j+1})) / (x_{j+3} - x_{j+1}),
 * где w1_j, w2_j — коэффициенты при производной из [DeBoorFixFunctionals] ([DerivFunctional.cD]).
 * Функционалы представлены как [ValueFunctional] (usesDerivative = false).
 *
 * Оператор P_xitilde — квазиинтерполянт, а не проектор (isProjector = false): при замене
 * производной разностью точная биортогональность в общем случае теряется. На span{1, rho, sigma}
 * погрешность имеет порядок O(h^2 |f''|); при кратных краевых узлах центральная разность на краю
 * становится односторонней, и в краевом слое порядок понижается. Конструкция собственная
 * (docs/REFERENCES.md, раздел 2).
 *
 * Краевые функционалы j = -2, n-1 — значения f(x_0), f(x_n), как в theta и xi.
 *
 * Вариант xitilde^{<0>} не реализуется: он опирался бы на f'' в левом конце носителя, где сплайн
 * и его первая производная обращаются в ноль, и замена второй производной разностью не
 * воспроизводит порождающую систему даже при h → 0.
 *
 * @property r порядок функционала из {1, 2}; по умолчанию r = 1.
 * @param ctx контекст численных вычислений; передаётся во вложенное семейство [DeBoorFixFunctionals].
 * @throws IllegalArgumentException если r вне {1, 2} либо вложенное семейство xi не строится.
 */
public class DiscreteDeBoorFixFunctionals(
    basis: MinimalSplineBasis,
    public val r: Int = 1,
    ctx: NumericsContext = NumericsContext.default(),
) : FunctionalFamily(basis, if (r == 1) "xitilde" else "xitilde<$r>", ctx) {
    init { require(r in 1..2) { "DiscreteDeBoorFix: r must be in {1,2}, got $r" } }
    override val isProjector: Boolean = false
    override val usesDerivative: Boolean = false
    // Контекст передаётся во вложенное семейство, чтобы `raw.ctx` совпадал с контекстом обёртки.
    private val raw = DeBoorFixFunctionals(basis, r, ctx)
    private val funcs: Array<ApproxFunctional> = Array(n + 2) { buildXiTilde(it - 2) }
    override fun chi(j: Int): ApproxFunctional = funcs[j + 2]

    private fun buildXiTilde(j: Int): ApproxFunctional {
        if (j == -2) return ValueFunctional(doubleArrayOf(grid.x(0)), doubleArrayOf(1.0))
        if (j == n - 1) return ValueFunctional(doubleArrayOf(grid.x(n)), doubleArrayOf(1.0))
        // Множитель при производной — тот же A^{<r>}_j, что и в исходном xi.
        val w = (raw.chi(j) as DerivFunctional).cD
        // Узел опоры и соседи для центральной разности вокруг него.
        val (node, left, right) = when (r) {
            2 -> Triple(grid.x(j + 2), grid.x(j + 1), grid.x(j + 3))
            else -> Triple(grid.x(j + 1), grid.x(j), grid.x(j + 2))
        }
        // Шаг разделённой разности — разность координат узлов; масштаб |right| + |left|.
        // Сама разность мала как h, поэтому абсолютный порог не отличал бы мелкую сетку от
        // совпавших узлов на отрезке, удалённом от нуля. См. KDoc splines.DEGENERACY_RELATIVE_EPS.
        val denom = right - left
        val denomScale = cancellationScale(right, left)
        require(isSignificant(denom, denomScale)) {
            "buildXiTilde(j=$j,r=$r): degenerate divided-difference span x=$right - x=$left = $denom, " +
                "scale=$denomScale (значимость потеряна: порог $DEGENERACY_RELATIVE_EPS)"
        }
        // f(node) + w (f(right) - f(left))/denom как комбинация значений.
        return ValueFunctional(
            doubleArrayOf(left, node, right),
            doubleArrayOf(-w / denom, 1.0, w / denom),
        )
    }
}

// ----------------------------------------------------------------------------
// mu — усредняющие функционалы
// Источник: Kulikov, Makarov, Construction of Approximation Functionals for
// Minimal Splines (Journal of Mathematical Sciences, 2022, vol. 262, no. 1,
// p. 84–98). См. docs/REFERENCES.md, раздел 2.
// ----------------------------------------------------------------------------

/**
 * Усредняющие функционалы mu_j(f) = a_j f(y_{j-1}) + b_j f(y_j) + c_j f(y_{j+1}).
 *
 * Узлы вспомогательной сетки: y_j = x_{j+1} + theta (x_{j+2} - x_{j+1}), по умолчанию
 * theta = 1/2 (середины интервалов). Коэффициенты определяются из условия точности на
 * span{1, rho, sigma}: решается система
 *   [1, 1, 1; rho(y_{j-1}), rho(y_j), rho(y_{j+1}); sigma(...)] (a, b, c)^T = a^N_j,
 * где a^N_j — тот же вектор аппроксимационного соотношения, который строит сам базис.
 *
 * Оператор P_mu — квазиинтерполянт, а не проектор (P_mu^2 != P_mu); точен на span{1, rho, sigma}.
 *
 * Контрольный частный случай: для полиномиальной системы B при theta = 1/2 на
 * равномерной сетке формула вырождается в -1/8 (f(y_{j-1}) - 10 f(y_j) + f(y_{j+1})).
 *
 * @property theta параметр размещения узлов вспомогательной сетки в (0, 1).
 * @param ctx контекст численных вычислений для СЛАУ 3×3.
 */
public class AveragingFunctionals(
    basis: MinimalSplineBasis,
    public val theta: Double = 0.5,
    ctx: NumericsContext = NumericsContext.default(),
) : FunctionalFamily(basis, "mu", ctx) {
    override val isProjector: Boolean = false
    override val usesDerivative: Boolean = false
    private val funcs: Array<ApproxFunctional> = Array(n + 2) { buildMu(it - 2) }
    override fun chi(j: Int): ApproxFunctional = funcs[j + 2]

    /** y_j по (net_Y): краевые y_{-2}=x_0, y_{n-1}=x_n; внутренние — x_{j+1}+theta(x_{j+2}-x_{j+1}). */
    private fun yNode(j: Int): Double = when (j) {
        -2 -> grid.x(0)
        n - 1 -> grid.x(n)
        else -> grid.x(j + 1) + theta * (grid.x(j + 2) - grid.x(j + 1))
    }

    private fun buildMu(j: Int): ApproxFunctional {
        if (j == -2) return ValueFunctional(doubleArrayOf(grid.x(0)), doubleArrayOf(1.0))
        if (j == n - 1) return ValueFunctional(doubleArrayOf(grid.x(n)), doubleArrayOf(1.0))
        val ym = yNode(j - 1); val y0 = yNode(j); val yp = yNode(j + 1)
        val ys = doubleArrayOf(ym, y0, yp)
        // Система (phi(y_{j-1}) | phi(y_j) | phi(y_{j+1})) mu = a_j записывается в локальных
        // координатах интервала (x_{j+1}, x_{j+2}) — среднего интервала носителя omega_j:
        // левое умножение обеих частей на T_{j+1} не меняет решения, а обусловленность матрицы
        // перестаёт зависеть от шага сетки и положения отрезка. Все три точки y_q лежат в
        // пределах двух шагов от x_{j+1}, где psi = O(1). Правая часть T_{j+1} a_j получается
        // той же формулой аппроксимационного соотношения, что и столбцы T_k M_k в базисе.
        val frame = basis.frame(j + 1)
        val cols = Array(3) { q -> frame.psi(ys[q]) }
        val matrix = DenseMatrix.build(3, 3) { r, c -> cols[c][r] }
        val coeff = LinearAlgebra.solve(matrix, basis.computeA(j, frame.psi, frame.psiD), ctx.backend)
        return ValueFunctional(ys, coeff)
    }
}

// ----------------------------------------------------------------------------
// lambda — трёхточечные функционалы
// Источник: Kulikov, Makarov (Journal of Mathematical Sciences, 2022, vol. 262,
// no. 1, p. 84–98). См. docs/REFERENCES.md, раздел 2.
// ----------------------------------------------------------------------------

/**
 * Трёхточечные функционалы lambda_j(f) по точкам x_{j+1}, x_{j+3/2}, x_{j+2},
 * где x_{j+3/2} = x_{j+1} + thetaHat (x_{j+2} - x_{j+1}).
 *
 * Реализованы через локальную аппроксимацию на отрезке I = [x_{j+1}, x_{j+2}]:
 * на нём активны ровно три сплайна omega_{j-1}, omega_j, omega_{j+1}; решается система
 * в трёх точках, и lambda_j(f) берётся как коэффициент при omega_j.
 *
 * Оператор P_lambda — квазиинтерполянт, а не проектор (P_lambda^2 != P_lambda); точен на span{1, rho, sigma}.
 *
 * Контрольный частный случай: для системы B при thetaHat = 1/2 формула вырождается
 * в -1/2 (f(x_{j+1}) - 4 f(x_{j+3/2}) + f(x_{j+2})).
 *
 * @property thetaHat параметр размещения средней точки в (0, 1); по умолчанию 1/2.
 * @param ctx контекст численных вычислений для СЛАУ 3×3.
 */
public class ThreePointFunctionals(
    basis: MinimalSplineBasis,
    public val thetaHat: Double = 0.5,
    ctx: NumericsContext = NumericsContext.default(),
) : FunctionalFamily(basis, "lambda", ctx) {
    override val isProjector: Boolean = false
    override val usesDerivative: Boolean = false
    private val funcs: Array<ApproxFunctional> = Array(n + 2) { buildLambda(it - 2) }
    override fun chi(j: Int): ApproxFunctional = funcs[j + 2]

    private fun buildLambda(j: Int): ApproxFunctional {
        if (j == -2) return ValueFunctional(doubleArrayOf(grid.x(0)), doubleArrayOf(1.0))
        if (j == n - 1) return ValueFunctional(doubleArrayOf(grid.x(n)), doubleArrayOf(1.0))
        val x1 = grid.x(j + 1); val x2 = grid.x(j + 2)
        val xMid = x1 + thetaHat * (x2 - x1)
        val points = doubleArrayOf(x1, xMid, x2)
        val active = intArrayOf(j - 1, j, j + 1) // активные на (x_{j+1},x_{j+2})
        // M[p][slot] = omega_active[slot](point_p); из M c = fvals следует lambda_j = c[1].
        // Коэффициенты coeff_p = (M^{-1})[1][p] — решение M^T coeff = e_1 (вторая строка обратной).
        val mTrans = DenseMatrix.build(3, 3) { i, p -> basis.omega(active[i], points[p]) }
        val e1 = doubleArrayOf(0.0, 1.0, 0.0)
        val coeff = LinearAlgebra.solve(mTrans, e1, ctx.backend)
        return ValueFunctional(points, coeff)
    }
}
