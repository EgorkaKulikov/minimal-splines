package splines.metrics

import kotlin.math.abs

import splines.Grid

// ============================================================================
// Метрики погрешности
// ============================================================================

/**
 * Кратность измельчения контрольной сетки в [errorEh] по умолчанию: `100 n` контрольных отрезков,
 * то есть `100 n + 1` контрольная точка.
 */
public const val DEFAULT_CONTROL_REFINEMENT: Int = 100

/**
 * Равномерная норма погрешности `E_h = max |exact(t) - eval(t)|` по `refinement * n + 1` точкам
 * равномерной контрольной сетки на отрезке `[a, b]` сетки [grid].
 *
 * Максимум по конечному набору точек не превышает точного максимума; достаточность значения
 * [refinement] по умолчанию проверяется пересчётом на более мелкой контрольной сетке.
 *
 * @param exact точная функция.
 * @param eval приближение (например, сплайн).
 * @param grid сетка, задающая отрезок `[a, b]` и число интервалов `n`.
 * @param refinement число контрольных отрезков на один интервал сетки.
 * @throws IllegalArgumentException если `grid.n < 1` или `refinement < 1`: число контрольных
 *   отрезков `m = refinement * n` служит делителем при вычислении контрольных точек.
 */
public fun errorEh(
    exact: (Double) -> Double,
    eval: (Double) -> Double,
    grid: Grid,
    refinement: Int = DEFAULT_CONTROL_REFINEMENT,
): Double {
    require(grid.n >= 1) { "errorEh: требуется grid.n >= 1 (иначе m = refinement*n = 0), получено n=${grid.n}" }
    require(refinement >= 1) { "errorEh: требуется refinement >= 1, получено $refinement" }
    val m = refinement * grid.n
    var e = 0.0
    for (i in 0..m) {
        // Произведение (b - a)·i/m при i = m может округлиться правее b на единицу
        // последнего разряда (например, b = 0.013, m = 130); сплайн вне [a, b] не определён,
        // поэтому контрольная точка ограничивается отрезком сетки. Для точек внутри
        // отрезка ограничение не меняет ни одного бита.
        val t = (grid.a + (grid.b - grid.a) * i / m).coerceIn(grid.a, grid.b)
        e = maxOf(e, abs(exact(t) - eval(t)))
    }
    return e
}
