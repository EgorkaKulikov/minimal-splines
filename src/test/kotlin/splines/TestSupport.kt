package splines

import splines.metrics.DEFAULT_CONTROL_REFINEMENT
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.exp
import kotlin.math.max
import kotlin.math.sin

/**
 * Максимум |exact − eval| по контрольной сетке `refinement·n + 1` точек, из которой исключён
 * краевой слой ширины `layer·h` с каждой стороны (h — максимальный шаг сетки).
 *
 * Нужен для семейства ξ̃ (`DiscreteDeBoorFixFunctionals`): у кратных краевых узлов центральная
 * разность становится односторонней, и квазипроектор не воспроизводит span φ в краевом слое
 * (дефект ≈ h²). Это свойство конструкции, а не дефект реализации, поэтому точность на span φ
 * и порядок 3 для ξ̃ проверяются во внутренней области.
 */
internal fun interiorErrorEh(
    exact: (Double) -> Double,
    eval: (Double) -> Double,
    grid: Grid,
    layer: Int = 3,
    refinement: Int = DEFAULT_CONTROL_REFINEMENT,
): Double {
    val m = refinement * grid.n
    val lo = grid.a + layer * grid.h
    val hi = grid.b - layer * grid.h
    var mx = 0.0
    for (i in 0..m) {
        val t = grid.a + (grid.b - grid.a) * i / m
        if (t < lo || t > hi) continue
        mx = max(mx, abs(exact(t) - eval(t)))
    }
    return mx
}

/** Тестовая функция f, f', f'' в виде тройки лямбд. */
internal data class TestFunction(
    val f: (Double) -> Double,
    val fD: (Double) -> Double,
    val fDD: (Double) -> Double,
)

/**
 * f(t) = exp(sin 3t), не лежащая в span φ ни одной из систем B, H, T;
 * f' = 3 cos 3t · f, f'' = (−9 sin 3t + 9 cos² 3t) · f.
 */
internal fun testFunction(): TestFunction = TestFunction(
    f = { t -> exp(sin(3.0 * t)) },
    fD = { t -> 3.0 * cos(3.0 * t) * exp(sin(3.0 * t)) },
    fDD = { t ->
        val s = sin(3.0 * t)
        val c = cos(3.0 * t)
        (-9.0 * s + 9.0 * c * c) * exp(s)
    },
)
