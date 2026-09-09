package splines

import numerics.Conditioning
import numerics.DenseMatrix
import numerics.NumericsContext

// ============================================================================
// 5. БАЗИС МИНИМАЛЬНЫХ СПЛАЙНОВ (M_k^{-1} phi(t) через numerical-core)
// ============================================================================

/**
 * Базис квадратичных минимальных сплайнов {omega_j}_{j=-2}^{n-1} на сетке с тройными краевыми узлами.
 * Значения на интервале (x_k, x_{k+1}) вычисляются как M_k^{-1} phi(t), где M_k = (a_{k-2}|a_{k-1}|a_k) —
 * матрица аппроксимационного соотношения; обратные матрицы вычисляются средствами numerical-core
 * в реализации BLAS/LAPACK из [ctx]. Вырожденность или недостоверность обращения приводит к исключению при построении.
 */
class MinimalSplineBasis(val sys: GeneratingSystem, val grid: Grid, ctx: NumericsContext = NumericsContext.default()) {
    val n = grid.n

    // a_j: на (x_{j+1},x_{j+2}) предел при тройном узле a_j = phi(x_{j+1}).
    private val aMin = -2
    private val aMax = n - 1
    private val aVec: Array<DoubleArray> = Array(aMax - aMin + 1) { k -> computeA(k + aMin) }

    private val invM: Array<DoubleArray> = Array(n) { k -> invertApproximationMatrix(k, ctx) }

    private fun a(j: Int): DoubleArray = aVec[j - aMin]

    /** Матрица аппроксимационного соотношения M_k = (a_{k-2}|a_{k-1}|a_k) на интервале (x_k, x_{k+1}). */
    internal fun approximationMatrix(k: Int): DenseMatrix =
        DenseMatrix.build(3, 3) { i, j -> a(k - 2 + j)[i] }

    private fun invertApproximationMatrix(k: Int, ctx: NumericsContext): DoubleArray {
        val m = approximationMatrix(k)
        val inv = Conditioning.inverse(m, ctx.backend)
            ?: throw IllegalArgumentException("Матрица аппроксимационного соотношения на интервале $k вырождена")
        val residual = Conditioning.inversionResidual(m, inv, ctx.backend)
        require(residual <= Conditioning.INVERSION_RESIDUAL_TOLERANCE) {
            "Матрица аппроксимационного соотношения на интервале $k плохо обусловлена: невязка обращения $residual"
        }
        return inv.data
    }

    /**
     * Вектор a_j аппроксимационного соотношения на (x_{j+1}, x_{j+2}), j = -2..n-1.
     *
     * `internal`, а не `private`: тот же вектор a^N_j нужен семейству усредняющих
     * функционалов mu (`AveragingFunctionals`), где раньше жила его дословная копия.
     */
    internal fun computeA(j: Int): DoubleArray {
        val xj1 = grid.x(j + 1)
        val phiJ1 = sys.phi(xj1)
        if (grid.isCoincident(j + 1)) return phiJ1 // тройной узел на краю: x_{j+1} = x_{j+2}
        val xj2 = grid.x(j + 2)
        val phiDJ1 = sys.phiD(xj1)
        val dJ2 = cross(sys.phi(xj2), sys.phiD(xj2))
        // Знаменатель — скалярное произведение, его масштаб задаёт сумма модулей
        // покомпонентных произведений (величина до взаимных сокращений). Абсолютный
        // порог здесь неприменим: denom ~ h, то есть зависит от масштаба отрезка.
        val denom = dot(dJ2, phiDJ1)
        val denomScale = dot3Scale(dJ2, phiDJ1)
        require(isSignificant(denom, denomScale)) {
            "computeA(j=$j): degenerate approximation relation, dot(dJ2, phiDJ1)=$denom, " +
                "scale=$denomScale (значимость потеряна: порог $DEGENERACY_RELATIVE_EPS)"
        }
        val coef = dot(dJ2, phiJ1) / denom
        return doubleArrayOf(
            phiJ1[0] - coef * phiDJ1[0],
            phiJ1[1] - coef * phiDJ1[1],
            phiJ1[2] - coef * phiDJ1[2],
        )
    }

    /**
     * Индекс сеточного интервала k с x_k <= t < x_{k+1} (для t=b возвращает n-1).
     *
     * ЕДИНСТВЕННАЯ ТОЧКА ПРОВЕРКИ ПРИНАДЛЕЖНОСТИ ОТРЕЗКУ для всего класса. Метод
     * `private`, но через него проходят ВСЕ пути, где выход за `[a,b]` вообще возможен:
     * [interval], [evalSpline], [evalSplineDeriv], [evalSplineDeriv2]. Семейство
     * `omega*` сюда с внешней точкой не попадает ПО ПОСТРОЕНИЮ: их предваряет отсечка
     * по носителю `x_j <= t <= x_{j+3}`, а носитель ЛЮБОГО базисного сплайна лежит
     * внутри отрезка — при j = -2 это `[x_{-2}, x_1] = [a, x_1]`, при j = n-1 это
     * `[x_{n-1}, x_{n+2}] = [x_{n-1}, b]` (краевые узлы трёхкратны). Поэтому проверка
     * здесь эквивалентна проверке в каждой публичной точке входа, но записана один раз.
     *
     * ЧТО БЫЛО РАНЬШЕ. Оба граничных условия ниже КЛАМПОВАЛИ внешнюю точку: `t < a`
     * давало интервал 0, `t > b` — интервал n-1. Формально это молчаливая ЭКСТРАПОЛЯЦИЯ:
     * [evalSpline] честно считал многочлен крайнего слоя в точке, где сплайн не
     * определён, и возвращал правдоподобное число вместо признака плохого входа.
     *
     * @throws IllegalArgumentException если t лежит строго вне [Grid.a], [Grid.b].
     */
    private fun intervalOf(t: Double): Int {
        // Условие записано ОТРИЦАНИЯМИ, а не как `t in grid.a..grid.b`, и это не стиль.
        // При t = NaN оба сравнения `<` и `>` ложны, поэтому оба отрицания истинны и NaN
        // проходит НАСКВОЗЬ — ровно как раньше (бинарный поиск ниже вернёт 1, phi(NaN)
        // даст NaN, и NaN распространится в результат). Это ТРЕБОВАНИЕ, а не побочный
        // эффект: `VolterraOperator.apply` сознательно пропускает NaN, чтобы плохой вход
        // проявился как NaN, а не как правдоподобное число (см.
        // `VolterraIntegrandCacheEquivalenceTest.nanArgument_bothPathsAgree`). Вариант
        // `t in grid.a..grid.b` на NaN дал бы false и превратил бы NaN в исключение,
        // сломав это соглашение.
        require(!(t < grid.a) && !(t > grid.b)) {
            val side = if (t < grid.a) "левее a на ${grid.a - t}" else "правее b на ${t - grid.b}"
            "MinimalSplineBasis: точка t=$t вне отрезка сетки [${grid.a}, ${grid.b}] ($side). " +
                "Сплайн определён ТОЛЬКО на отрезке сетки: вне его значение не задано ни " +
                "аппроксимационным соотношением, ни носителями базиса. Прежнее поведение " +
                "МОЛЧА экстраполировало многочлен крайнего слоя и возвращало правдоподобное " +
                "число; теперь такой вызов — ошибка. Если точка получена шаблоном конечной " +
                "разности или квадратурой по [a,t] при t>b, ограничьте аргумент отрезком явно."
        }
        // Бинарный поиск (breakpoints x_0..x_n возрастают): наибольший k in [0,n-1]
        // с x_k <= t, с клампами на концах. Семантика идентична линейному поиску:
        // t < x_1 -> 0; t >= x_{n-1} -> n-1; иначе x_lo <= t < x_{lo+1}.
        if (t < grid.x(1)) return 0
        if (t >= grid.x(n - 1)) return n - 1
        var lo = 1
        var hi = n - 1
        while (hi - lo > 1) {
            val mid = (lo + hi) / 2
            if (grid.x(mid) <= t) lo = mid else hi = mid
        }
        return lo
    }

    /**
     * Индекс сеточного интервала, содержащего t (публичный доступ).
     *
     * @throws IllegalArgumentException если t строго вне отрезка сетки (см. [intervalOf]).
     */
    fun interval(t: Double): Int = intervalOf(t)

    /** Три активных значения omega_{k-2},omega_{k-1},omega_k в точке t (одно M_k^{-1} phi(t)). */
    fun activeOmega(k: Int, t: Double): DoubleArray {
        val inv = invM[k]
        val p = sys.phi(t)
        return doubleArrayOf(
            inv[0] * p[0] + inv[3] * p[1] + inv[6] * p[2],
            inv[1] * p[0] + inv[4] * p[1] + inv[7] * p[2],
            inv[2] * p[0] + inv[5] * p[1] + inv[8] * p[2],
        )
    }

    /** Значение omega_j(t), j из [-2, n-1]. Носитель [x_j, x_{j+3}]; вне него 0. */
    fun omega(j: Int, t: Double): Double {
        if (t < grid.x(j) || t > grid.x(j + 3)) return 0.0
        val k = intervalOf(t)
        val slot = j - (k - 2)
        if (slot < 0 || slot > 2) return 0.0
        val inv = invM[k]
        val phiT = sys.phi(t)
        return inv[slot] * phiT[0] + inv[slot + 3] * phiT[1] + inv[slot + 6] * phiT[2]
    }

    /** Производная omega_j'(t) (phi заменяется на phi'). Нужна для xi-функционалов. */
    fun omegaDeriv(j: Int, t: Double): Double {
        if (t < grid.x(j) || t > grid.x(j + 3)) return 0.0
        val k = intervalOf(t)
        val slot = j - (k - 2)
        if (slot < 0 || slot > 2) return 0.0
        val inv = invM[k]
        val phiDT = sys.phiD(t)
        return inv[slot] * phiDT[0] + inv[slot + 3] * phiDT[1] + inv[slot + 6] * phiDT[2]
    }

    /**
     * Вторая производная omega_j''(t) (phi заменяется на phi''). Нужна для xi^<0>
     * (де Бура--Фикса r=0). Кусочно-постоянна по слоям; в узлах сетки omega_j'' терпит
     * разрыв (omega_j in C^1 \ C^2), поэтому значение в узле берётся по правому куску.
     */
    fun omegaDeriv2(j: Int, t: Double): Double {
        if (t < grid.x(j) || t > grid.x(j + 3)) return 0.0
        val k = intervalOf(t)
        val slot = j - (k - 2)
        if (slot < 0 || slot > 2) return 0.0
        val inv = invM[k]
        val phiDDT = sys.phiDD(t)
        return inv[slot] * phiDDT[0] + inv[slot + 3] * phiDDT[1] + inv[slot + 6] * phiDDT[2]
    }

    /**
     * Значение сплайна u_h(t) = sum_j c_j omega_j(t), c размера n+2.
     *
     * @throws IllegalArgumentException если t строго вне отрезка сетки (см. [intervalOf]):
     *   вне отрезка сплайн не определён, а прежнее поведение молча экстраполировало.
     */
    fun evalSpline(c: DoubleArray, t: Double): Double {
        val k = intervalOf(t)
        val w = activeOmega(k, t)
        return c[k] * w[0] + c[k + 1] * w[1] + c[k + 2] * w[2] // индексы k-2,k-1,k -> +2
    }

    /**
     * Значение производной сплайна u_h'(t).
     *
     * @throws IllegalArgumentException если t строго вне отрезка сетки (см. [intervalOf]).
     */
    fun evalSplineDeriv(c: DoubleArray, t: Double): Double {
        val k = intervalOf(t)
        val inv = invM[k]
        val p = sys.phiD(t)
        val w0 = inv[0] * p[0] + inv[3] * p[1] + inv[6] * p[2]
        val w1 = inv[1] * p[0] + inv[4] * p[1] + inv[7] * p[2]
        val w2 = inv[2] * p[0] + inv[5] * p[1] + inv[8] * p[2]
        return c[k] * w0 + c[k + 1] * w1 + c[k + 2] * w2
    }

    /**
     * Значение второй производной сплайна u_h''(t) (для xi^{<0>}-идемпотентности).
     *
     * @throws IllegalArgumentException если t строго вне отрезка сетки (см. [intervalOf]).
     */
    fun evalSplineDeriv2(c: DoubleArray, t: Double): Double {
        val k = intervalOf(t)
        val inv = invM[k]
        val p = sys.phiDD(t)
        val w0 = inv[0] * p[0] + inv[3] * p[1] + inv[6] * p[2]
        val w1 = inv[1] * p[0] + inv[4] * p[1] + inv[7] * p[2]
        val w2 = inv[2] * p[0] + inv[5] * p[1] + inv[8] * p[2]
        return c[k] * w0 + c[k + 1] * w1 + c[k + 2] * w2
    }
}

// ----------------------------------------------------------------------------
// 5.1 ЭТАЛОННЫЕ ФОРМУЛЫ базисов (для health-checks)
// ----------------------------------------------------------------------------

/** Эталонные явные формулы для сверки с общим базисом (health-check). */
object ReferenceSplines {
    /** Классический квадратичный B-сплайн omega^B_j(t). */
    fun omegaB(grid: Grid, j: Int, t: Double): Double {
        val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2); val xj3 = grid.x(j + 3)
        if (t < xj || t > xj3) return 0.0
        return when {
            t < xj1 -> (t - xj) * (t - xj) / ((xj1 - xj) * (xj2 - xj))
            t < xj2 -> (1.0 / (xj1 - xj)) * (
                (t - xj) * (t - xj) / (xj2 - xj)
                    - (t - xj1) * (t - xj1) * (xj3 - xj) / ((xj2 - xj1) * (xj3 - xj1))
                )
            else -> (t - xj3) * (t - xj3) / ((xj3 - xj1) * (xj3 - xj2))
        }
    }

    /** Производная эталонного B-сплайна. */
    fun omegaBDeriv(grid: Grid, j: Int, t: Double): Double {
        val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2); val xj3 = grid.x(j + 3)
        if (t < xj || t > xj3) return 0.0
        return when {
            t < xj1 -> 2.0 * (t - xj) / ((xj1 - xj) * (xj2 - xj))
            t < xj2 -> (1.0 / (xj1 - xj)) * (
                2.0 * (t - xj) / (xj2 - xj)
                    - 2.0 * (t - xj1) * (xj3 - xj) / ((xj2 - xj1) * (xj3 - xj1))
                )
            else -> 2.0 * (t - xj3) / ((xj3 - xj1) * (xj3 - xj2))
        }
    }

    /** Гиперболический минимальный сплайн omega^H_j(t). */
    fun omegaH(grid: Grid, j: Int, t: Double): Double {
        val xj = grid.x(j); val xj1 = grid.x(j + 1); val xj2 = grid.x(j + 2); val xj3 = grid.x(j + 3)
        if (t < xj || t > xj3) return 0.0
        fun sh(x: Double) = Math.sinh(x)
        fun ch(x: Double) = Math.cosh(x)
        val cosCenter = ch((xj2 - xj1) / 2.0)
        return when {
            t < xj1 -> cosCenter * sh((t - xj) / 2.0) * sh((t - xj) / 2.0) /
                (sh((xj1 - xj) / 2.0) * sh((xj2 - xj) / 2.0))
            t < xj2 -> (cosCenter / sh((xj1 - xj) / 2.0)) * (
                sh((t - xj) / 2.0) * sh((t - xj) / 2.0) / sh((xj2 - xj) / 2.0)
                    - sh((xj3 - xj) / 2.0) * sh((t - xj1) / 2.0) * sh((t - xj1) / 2.0) /
                    (sh((xj3 - xj1) / 2.0) * sh((xj2 - xj1) / 2.0))
                )
            else -> cosCenter * sh((xj3 - t) / 2.0) * sh((xj3 - t) / 2.0) /
                (sh((xj3 - xj1) / 2.0) * sh((xj3 - xj2) / 2.0))
        }
    }
}

/** Узлы x_j..x_{j+3} различны (нет слияния кратных узлов). */
fun nonDegenerate(grid: Grid, j: Int): Boolean =
    grid.x(j) < grid.x(j + 1) && grid.x(j + 1) < grid.x(j + 2) && grid.x(j + 2) < grid.x(j + 3)

/** Векторное произведение u x v в R^3 — локальная арифметика построения a_j. */
private fun cross(u: DoubleArray, v: DoubleArray): DoubleArray = doubleArrayOf(
    u[1] * v[2] - u[2] * v[1],
    u[2] * v[0] - u[0] * v[2],
    u[0] * v[1] - u[1] * v[0],
)

/** Скалярное произведение в R^3. */
private fun dot(u: DoubleArray, v: DoubleArray): Double = u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
