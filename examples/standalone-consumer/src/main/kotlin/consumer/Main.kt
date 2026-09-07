package consumer

import numerics.orders
import numerics.reliableOrders
import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import splines.functionals.AveragingFunctionals
import splines.functionals.FunctionalFamily
import splines.functionals.ProjFunctionals
import splines.metrics.errorEh
import kotlin.math.cos
import kotlin.math.exp
import kotlin.math.sin

/**
 * Сторонний потребитель библиотек `numerical-core` и `minimal-splines`.
 *
 * Численный метод, не связанный с интегральными уравнениями: квазиинтерполяция гладкой
 * функции f(t) = exp(sin 3t) квадратичными минимальными сплайнами и численное
 * дифференцирование через производную сплайна. Импортируются ТОЛЬКО пакеты `numerics.*`
 * и `splines.*` — так фиксируется контракт: для работы со сплайнами репозиторий
 * `integral-equations` не нужен.
 */
fun f(t: Double): Double = exp(sin(3.0 * t))
fun fD(t: Double): Double = 3.0 * cos(3.0 * t) * exp(sin(3.0 * t))

private val ns = listOf(8, 16, 32, 64)

private fun table(title: String, family: (MinimalSplineBasis) -> FunctionalFamily) {
    val errs = ArrayList<Double>()
    val derivErrs = ArrayList<Double>()
    for (n in ns) {
        val grid = Grid.uniform(n, 0.0, 1.0)
        val basis = MinimalSplineBasis(GeneratingSystem.B, grid)
        val funcs = family(basis)
        // Коэффициенты (квази)проекции P_chi f = sum_j chi_j(f) omega_j.
        val c = funcs.projectorCoeffs(::f)
        errs += errorEh(::f, { t -> basis.evalSpline(c, t) }, grid)
        derivErrs += errorEh(::fD, { t -> basis.evalSplineDeriv(c, t) }, grid)
    }
    val p = orders(errs)
    val pReliable = reliableOrders(errs)
    val pD = orders(derivErrs)
    println(title)
    println("  n     E_h          p_h     p_h(reliable)   E_h(f')      p_h(f')")
    for (i in ns.indices) {
        println(
            "  %-4d  %-11.3e  %-6s  %-14s  %-11.3e  %s".format(
                ns[i], errs[i], fmt(p[i]), pReliable[i]?.let { fmt(it) } ?: "---", derivErrs[i], fmt(pD[i]),
            ),
        )
    }
    println()
}

private fun fmt(v: Double): String = if (v.isNaN()) "---" else "%.2f".format(v)

fun main() {
    println("Квазиинтерполяция f(t) = exp(sin 3t) на [0, 1] квадратичными минимальными сплайнами (система B)")
    println("E_h — равномерная норма ошибки на контрольной сетке, p_h = log2(E_h / E_{h/2})")
    println()
    table("theta — проекционные функционалы (ProjFunctionals)") { ProjFunctionals(it) }
    table("mu — усредняющие функционалы (AveragingFunctionals)") { AveragingFunctionals(it) }

    // Численное дифференцирование напрямую через производные базисных сплайнов:
    // f'(t) ~ sum_j chi_j(f) omega_j'(t). Проверяем разбиение единицы и её производную.
    val grid = Grid.uniform(16, 0.0, 1.0)
    val basis = MinimalSplineBasis(GeneratingSystem.B, grid)
    val t = 0.37
    var sum = 0.0
    var sumD = 0.0
    for (j in -2..grid.n - 1) {
        sum += basis.omega(j, t)
        sumD += basis.omegaDeriv(j, t)
    }
    println("Разбиение единицы в t=$t: sum_j omega_j = %.15f, sum_j omega_j' = %.2e".format(sum, sumD))
}
