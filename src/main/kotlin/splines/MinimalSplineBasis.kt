package splines

import numerics.Conditioning
import numerics.DenseMatrix
import numerics.NumericsContext

// ============================================================================
// Базис минимальных сплайнов: omega(t) = M_k^{-1} phi(t) в локальных координатах интервала
// ============================================================================

/**
 * Базис квадратичных минимальных сплайнов {omega_j}_{j=-2}^{n-1} на сетке с тройными краевыми узлами.
 *
 * На интервале (x_k, x_{k+1}) значения задаются аппроксимационным соотношением
 * sum_j a_j omega_j(t) = phi(t), то есть omega(t) = M_k^{-1} phi(t), где M_k = (a_{k-2}|a_{k-1}|a_k).
 * Вычисления ведутся в локальном представлении порождающей системы psi_k = T_k phi
 * (см. [GeneratingSystem.localFrame]): левое умножение на невырожденную матрицу T_k даёт
 * эквивалентную систему (T_k M_k) omega(t) = psi_k(t) с тем же решением, поэтому базисные
 * функции от выбора представления не зависят. Столбцы T_k a_j вычисляются формулой [computeA],
 * записанной для psi_k вместо phi: вектор a_j ковариантен относительно замены phi -> T phi
 * (a_j^{T phi} = T a_j^{phi}), так как коэффициент при phi'(x_{j+1}) есть отношение двух скалярных
 * произведений с одним и тем же вектором нормали и при линейной замене не меняется.
 * Компоненты psi_k на интервале имеют порядок единицы, и число обусловленности T_k M_k не зависит
 * ни от числа интервалов, ни от положения и масштаба отрезка.
 *
 * Обратные матрицы вычисляются реализацией BLAS/LAPACK из контекста ctx; достоверность
 * обращения контролируется оценкой числа обусловленности T_k M_k (см. [MAX_CONDITION]).
 *
 * @property sys порождающая система phi = (1, rho, sigma).
 * @property grid сетка с кратными краевыми узлами.
 * @param ctx контекст численных вычислений (реализация BLAS/LAPACK).
 * @throws IllegalArgumentException если на некотором интервале порождающая система переполняется
 *   в double либо матрица аппроксимационного соотношения вырождена или её число обусловленности
 *   превышает [MAX_CONDITION].
 */
public class MinimalSplineBasis(public val sys: GeneratingSystem, public val grid: Grid, ctx: NumericsContext = NumericsContext.default()) {
    /** Константы критерия вырожденности. */
    public companion object {
        /**
         * Наибольшее допустимое число обусловленности матрицы аппроксимационного соотношения
         * в локальных координатах интервала; при большем значении в обратной матрице сохраняется
         * менее восьми значащих цифр, и базис считается вырожденным на данном интервале.
         */
        public const val MAX_CONDITION: Double = 1e8
    }

    /** Число интервалов сетки; базис состоит из n + 2 функций omega_{-2}, ..., omega_{n-1}. */
    public val n: Int = grid.n

    /** Локальные представления порождающей системы на интервалах (x_k, x_{k+1}), k = 0..n-1. */
    private val frames: Array<LocalFrame> = Array(n) { k -> sys.localFrame(grid.x(k), grid.x(k + 1) - grid.x(k)) }

    /** Локальное представление psi_k = T_k phi для интервала (x_k, x_{k+1}); используется семействами функционалов. */
    internal fun frame(k: Int): LocalFrame = frames[k]

    /** Обратные матрицы (T_k M_k)^{-1} по столбцам: элемент (slot, p) хранится в data[slot + 3 p]. */
    private val invM: Array<DoubleArray> = Array(n) { k -> invertApproximationMatrix(k, ctx) }

    /** Матрица аппроксимационного соотношения M_k = (a_{k-2}|a_{k-1}|a_k) в глобальных координатах phi. */
    internal fun approximationMatrix(k: Int): DenseMatrix {
        val cols = Array(3) { j -> computeA(k - 2 + j) }
        return DenseMatrix.build(3, 3) { i, j -> cols[j][i] }
    }

    /** Матрица T_k M_k = (T_k a_{k-2}|T_k a_{k-1}|T_k a_k) в локальных координатах интервала k; именно она обращается. */
    internal fun localApproximationMatrix(k: Int): DenseMatrix {
        val frame = frames[k]
        requireFiniteFrame(k, frame)
        val cols = Array(3) { j -> computeA(k - 2 + j, frame.psi, frame.psiD) }
        return DenseMatrix.build(3, 3) { i, j -> cols[j][i] }
    }

    /**
     * Проверка представимости порождающей системы в double на интервале k: столбцы a_{k-2}, a_{k-1}, a_k
     * строятся по значениям psi_k, psi_k' в узлах x_{k-1}, ..., x_{k+2}, и все они должны быть финитны.
     * Для систем H и T с локальным масштабом l = min(h, 1) значения sinh(u/l) переполняются уже при
     * u/l > 710, то есть при шаге h > 355 для H; без этой проверки переполнение проявлялось бы как
     * NaN в сообщении о вырожденности аппроксимационного соотношения.
     */
    private fun requireFiniteFrame(k: Int, frame: LocalFrame) {
        for (m in k - 1..k + 2) {
            val x = grid.x(m)
            require(allFinite(frame.psi(x)) && allFinite(frame.psiD(x))) {
                "Порождающая система ${sys.name} переполняется на интервале $k = [${grid.x(k)}, ${grid.x(k + 1)}]: " +
                    "значения psi_k или psi_k' в узле x_$m = $x нефинитны; уменьшите шаг сетки или длину отрезка"
            }
        }
    }

    private fun invertApproximationMatrix(k: Int, ctx: NumericsContext): DoubleArray {
        val m = localApproximationMatrix(k)
        val cond = Conditioning.conditionEstimate(m, ctx).valueOrNull()
            ?: throw IllegalArgumentException(
                "Матрица аппроксимационного соотношения на интервале $k численно вырождена: " +
                    "оценка числа обусловленности недостоверна",
            )
        if (cond > MAX_CONDITION) {
            throw IllegalArgumentException(
                "Матрица аппроксимационного соотношения на интервале $k вырождена: " +
                    "число обусловленности $cond превышает MAX_CONDITION = $MAX_CONDITION",
            )
        }
        val inv = Conditioning.inverse(m, ctx.backend)
            ?: throw IllegalArgumentException(
                "Матрица аппроксимационного соотношения на интервале $k численно вырождена",
            )
        return inv.data
    }

    /**
     * Вектор a_j аппроксимационного соотношения на (x_{j+1}, x_{j+2}), j = -2..n-1,
     * в глобальных координатах phi.
     *
     * `internal`, а не `private`: тот же вектор a^N_j нужен семейству усредняющих
     * функционалов mu (`AveragingFunctionals`). В построении базиса не участвует: столбцы T_k M_k
     * вычисляются через локальное представление [frame]. Для систем H и T на отрезках с |t| > 700
     * глобальные значения sinh, cosh переполняются, и метод неприменим (см. [requireFiniteFrame]).
     */
    internal fun computeA(j: Int): DoubleArray = computeA(j, sys::phi, sys::phiD)

    /**
     * Вектор a_j для порождающей вектор-функции [phi] с производной [phiD]: направляющий вектор
     * пересечения плоскостей span{phi(x_{j+1}), phi'(x_{j+1})} и span{phi(x_{j+2}), phi'(x_{j+2})},
     * нормированный условием a_j = phi(x_{j+1}) - coef phi'(x_{j+1}). Коэффициент coef — отношение
     * скалярных произведений с общей нормалью phi(x_{j+2}) × phi'(x_{j+2}), поэтому при замене
     * phi -> T phi вектор переходит в T a_j; это позволяет вычислять столбцы T_k M_k той же формулой.
     */
    internal fun computeA(j: Int, phi: (Double) -> DoubleArray, phiD: (Double) -> DoubleArray): DoubleArray {
        val xj1 = grid.x(j + 1)
        val phiJ1 = phi(xj1)
        if (grid.isCoincident(j + 1)) return phiJ1 // тройной узел на краю: x_{j+1} = x_{j+2}
        val xj2 = grid.x(j + 2)
        val phiDJ1 = phiD(xj1)
        val dJ2 = cross(phi(xj2), phiD(xj2))
        // Знаменатель — скалярное произведение, его масштаб задаёт сумма модулей
        // покомпонентных произведений (величина до взаимных сокращений). Абсолютный
        // порог здесь неприменим: denom ~ h, то есть зависит от масштаба отрезка.
        val denom = dot(dJ2, phiDJ1)
        val denomScale = dot3Scale(dJ2, phiDJ1)
        require(denomScale.isFinite()) {
            "computeA(j=$j): произведения значений phi и phi' в узлах x_{j+1}=$xj1, x_{j+2}=$xj2 " +
                "переполняются (масштаб $denomScale); уменьшите шаг сетки или длину отрезка"
        }
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
     * Индекс сеточного интервала k с x_k <= t < x_{k+1} (для t = b возвращается n-1).
     *
     * Единственная точка проверки принадлежности точки отрезку: через неё проходят [interval],
     * [evalSpline], [evalSplineDeriv], [evalSplineDeriv2]. Методы omega* с внешней точкой сюда не
     * обращаются: их предваряет отсечка по носителю [x_j, x_{j+3}], лежащему внутри отрезка.
     *
     * @throws IllegalArgumentException если t лежит строго вне [Grid.a], [Grid.b].
     */
    private fun intervalOf(t: Double): Int {
        // Условие записано через отрицания, а не как `t in grid.a..grid.b`: при t = NaN оба
        // сравнения ложны, и NaN проходит далее, распространяясь в результат (phi(NaN) = NaN);
        // запись через диапазон превратила бы NaN в исключение.
        require(!(t < grid.a) && !(t > grid.b)) {
            val side = if (t < grid.a) "левее a на ${grid.a - t}" else "правее b на ${t - grid.b}"
            "MinimalSplineBasis: точка t=$t вне отрезка сетки [${grid.a}, ${grid.b}] ($side). " +
                "Сплайн определён только на отрезке сетки, и вне его значение не экстраполируется. " +
                "Если точка получена шаблоном конечной разности или квадратурой при t > b, " +
                "аргумент следует ограничить отрезком."
        }
        // Бинарный поиск (breakpoints x_0..x_n возрастают): наибольший k in [0,n-1]
        // с x_k <= t, с ограничением на концах. Результат совпадает с линейным поиском:
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
    public fun interval(t: Double): Int = intervalOf(t)

    /** Проверка индекса базисного сплайна j из [-2, n-1]. */
    private fun requireIndex(j: Int) {
        require(j in -2..n - 1) { "Индекс сплайна j должен лежать в [-2, ${n - 1}], получено $j" }
    }

    /** Проверка длины вектора коэффициентов: по одному на базисный сплайн omega_{-2}, ..., omega_{n-1}. */
    private fun requireCoefficients(c: DoubleArray) {
        require(c.size == n + 2) { "Вектор коэффициентов должен иметь длину n + 2 = ${n + 2}, получено ${c.size}" }
    }

    /**
     * Три активных значения omega_{k-2},omega_{k-1},omega_k в точке t (одно (T_k M_k)^{-1} psi_k(t)).
     *
     * @throws IllegalArgumentException если k вне [0, n-1].
     */
    public fun activeOmega(k: Int, t: Double): DoubleArray {
        require(k in 0..n - 1) { "Индекс интервала k должен лежать в [0, ${n - 1}], получено $k" }
        val inv = invM[k]
        val p = frames[k].psi(t)
        return doubleArrayOf(
            inv[0] * p[0] + inv[3] * p[1] + inv[6] * p[2],
            inv[1] * p[0] + inv[4] * p[1] + inv[7] * p[2],
            inv[2] * p[0] + inv[5] * p[1] + inv[8] * p[2],
        )
    }

    /**
     * Значение omega_j(t), j из [-2, n-1]. Носитель [x_j, x_{j+3}]; вне него 0.
     *
     * @throws IllegalArgumentException если j вне [-2, n-1].
     */
    public fun omega(j: Int, t: Double): Double {
        requireIndex(j)
        if (t < grid.x(j) || t > grid.x(j + 3)) return 0.0
        val k = intervalOf(t)
        val slot = j - (k - 2)
        if (slot < 0 || slot > 2) return 0.0
        val inv = invM[k]
        val p = frames[k].psi(t)
        return inv[slot] * p[0] + inv[slot + 3] * p[1] + inv[slot + 6] * p[2]
    }

    /**
     * Производная omega_j'(t) (psi_k заменяется на psi_k'). Нужна для xi-функционалов.
     *
     * @throws IllegalArgumentException если j вне [-2, n-1].
     */
    public fun omegaDeriv(j: Int, t: Double): Double {
        requireIndex(j)
        if (t < grid.x(j) || t > grid.x(j + 3)) return 0.0
        val k = intervalOf(t)
        val slot = j - (k - 2)
        if (slot < 0 || slot > 2) return 0.0
        val inv = invM[k]
        val p = frames[k].psiD(t)
        return inv[slot] * p[0] + inv[slot + 3] * p[1] + inv[slot + 6] * p[2]
    }

    /**
     * Вторая производная omega_j''(t) (psi_k заменяется на psi_k''). Нужна для xi^<0>
     * (де Бура--Фикса r=0). Кусочно-постоянна по слоям; в узлах сетки omega_j'' терпит
     * разрыв (omega_j in C^1 \ C^2), поэтому значение в узле берётся по правому куску.
     *
     * @throws IllegalArgumentException если j вне [-2, n-1].
     */
    public fun omegaDeriv2(j: Int, t: Double): Double {
        requireIndex(j)
        if (t < grid.x(j) || t > grid.x(j + 3)) return 0.0
        val k = intervalOf(t)
        val slot = j - (k - 2)
        if (slot < 0 || slot > 2) return 0.0
        val inv = invM[k]
        val p = frames[k].psiDD(t)
        return inv[slot] * p[0] + inv[slot + 3] * p[1] + inv[slot + 6] * p[2]
    }

    /**
     * Значение сплайна u_h(t) = sum_j c_j omega_j(t), c размера n+2.
     *
     * @throws IllegalArgumentException если c.size != n + 2 или t строго вне отрезка сетки (см. [intervalOf]).
     */
    public fun evalSpline(c: DoubleArray, t: Double): Double {
        requireCoefficients(c)
        val k = intervalOf(t)
        val w = activeOmega(k, t)
        return c[k] * w[0] + c[k + 1] * w[1] + c[k + 2] * w[2] // индексы k-2,k-1,k -> +2
    }

    /**
     * Значение производной сплайна u_h'(t).
     *
     * @throws IllegalArgumentException если c.size != n + 2 или t строго вне отрезка сетки (см. [intervalOf]).
     */
    public fun evalSplineDeriv(c: DoubleArray, t: Double): Double {
        requireCoefficients(c)
        val k = intervalOf(t)
        val inv = invM[k]
        val p = frames[k].psiD(t)
        val w0 = inv[0] * p[0] + inv[3] * p[1] + inv[6] * p[2]
        val w1 = inv[1] * p[0] + inv[4] * p[1] + inv[7] * p[2]
        val w2 = inv[2] * p[0] + inv[5] * p[1] + inv[8] * p[2]
        return c[k] * w0 + c[k + 1] * w1 + c[k + 2] * w2
    }

    /**
     * Значение второй производной сплайна u_h''(t) (для xi^{<0>}-идемпотентности).
     *
     * @throws IllegalArgumentException если c.size != n + 2 или t строго вне отрезка сетки (см. [intervalOf]).
     */
    public fun evalSplineDeriv2(c: DoubleArray, t: Double): Double {
        requireCoefficients(c)
        val k = intervalOf(t)
        val inv = invM[k]
        val p = frames[k].psiDD(t)
        val w0 = inv[0] * p[0] + inv[3] * p[1] + inv[6] * p[2]
        val w1 = inv[1] * p[0] + inv[4] * p[1] + inv[7] * p[2]
        val w2 = inv[2] * p[0] + inv[5] * p[1] + inv[8] * p[2]
        return c[k] * w0 + c[k + 1] * w1 + c[k + 2] * w2
    }
}

/** Узлы x_j..x_{j+3} различны (нет слияния кратных узлов). */
public fun nonDegenerate(grid: Grid, j: Int): Boolean =
    grid.x(j) < grid.x(j + 1) && grid.x(j + 1) < grid.x(j + 2) && grid.x(j + 2) < grid.x(j + 3)

/** Векторное произведение u x v в R^3 — локальная арифметика построения a_j. */
private fun cross(u: DoubleArray, v: DoubleArray): DoubleArray = doubleArrayOf(
    u[1] * v[2] - u[2] * v[1],
    u[2] * v[0] - u[0] * v[2],
    u[0] * v[1] - u[1] * v[0],
)

/** Скалярное произведение в R^3. */
private fun dot(u: DoubleArray, v: DoubleArray): Double = u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

/** Все компоненты вектора финитны (нет переполнения и NaN). */
private fun allFinite(u: DoubleArray): Boolean = u.all { it.isFinite() }
