package splines.bench

import numerics.backend.Backends
import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import splines.functionals.AveragingFunctionals
import splines.functionals.DeBoorFixFunctionals
import splines.functionals.ProjFunctionals
import java.time.OffsetDateTime
import java.time.format.DateTimeFormatter
import java.util.Random
import kotlin.math.cos
import kotlin.math.exp
import kotlin.math.sin

// Измерения производительности публичного API: построение базиса, три семейства
// функционалов и вычисление сплайна. Не тест — запускается ./gradlew benchmark.

private fun f(t: Double) = exp(sin(3 * t))
private fun fD(t: Double) = 3 * cos(3 * t) * exp(sin(3 * t))
private fun fDD(t: Double) = (-9 * sin(3 * t) + 9 * cos(3 * t) * cos(3 * t)) * exp(sin(3 * t))

/** Медиана из [reps] замеров (мс) после [warm] прогревов. */
fun median(warm: Int, reps: Int, block: () -> Unit): Double {
    repeat(warm) { block() }
    val t = DoubleArray(reps) {
        val s = System.nanoTime()
        block()
        (System.nanoTime() - s) / 1e6
    }
    t.sort()
    return t[reps / 2]
}

/** Один замер: если он дольше 2 с — медиана из 3 после 1 прогрева, иначе из 5 после 3. */
private fun timed(block: () -> Unit): Double {
    val s = System.nanoTime()
    block()
    val first = (System.nanoTime() - s) / 1e6
    return if (first > 2000.0) median(1, 3, block) else median(3, 5, block)
}

private fun fmt(ms: Double) = String.format("%.2f", ms)

fun main(args: Array<String>) {
    val sizes = (if (args.isEmpty()) listOf("100", "1000", "10000") else args.toList()).map { it.toInt() }
    println("backend: ${Backends.default().name}")
    println("processors: ${Runtime.getRuntime().availableProcessors()}")
    println("jvm: ${System.getProperty("java.vm.name")} ${System.getProperty("java.version")}")
    println("date: ${OffsetDateTime.now().format(DateTimeFormatter.ISO_OFFSET_DATE_TIME)}")
    println()
    println("| система | n | basis | theta | xi1 | mu | eval×10⁴ |")
    println("|---|---|---|---|---|---|---|")
    val systems = listOf(GeneratingSystem.B, GeneratingSystem.H)
    val evalCount = 10_000
    for (n in sizes) {
        for (sys in systems) {
            val grid = Grid.uniform(n, 0.0, 1.0)
            // Побочный результат нужен, чтобы JIT не выбросил вычисление.
            var sink = 0.0
            // Базис может не строиться по критерию обусловленности (cond(M_k) > MAX_CONDITION):
            // такая строка выводится с прочерками и причиной, остальные размеры измеряются.
            val basis = try {
                MinimalSplineBasis(sys, grid)
            } catch (e: IllegalArgumentException) {
                println("| ${sys.name} | $n | — | — | — | — | — | ${e.message} |")
                continue
            }
            val tBasis = timed { sink += MinimalSplineBasis(sys, grid).grid.n }
            val tTheta = timed { sink += ProjFunctionals(basis).projectorCoeffs(::f, ::fD, ::fDD)[0] }
            val tXi1 = timed { sink += DeBoorFixFunctionals(basis, 1).projectorCoeffs(::f, ::fD, ::fDD)[0] }
            val tMu = timed { sink += AveragingFunctionals(basis).projectorCoeffs(::f, ::fD, ::fDD)[0] }
            val c = ProjFunctionals(basis).projectorCoeffs(::f, ::fD, ::fDD)
            val ts = Random(42).let { rnd -> DoubleArray(evalCount) { rnd.nextDouble() } }
            val tEval = timed { for (t in ts) sink += basis.evalSpline(c, t) }
            if (sink.isNaN()) println("nan")
            println("| ${sys.name} | $n | ${fmt(tBasis)} | ${fmt(tTheta)} | ${fmt(tXi1)} | ${fmt(tMu)} | ${fmt(tEval)} |")
        }
    }
}
