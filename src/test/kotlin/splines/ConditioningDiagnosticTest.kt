package splines

import numerics.Conditioning
import numerics.backend.Backends
import org.junit.jupiter.api.Tag
import java.io.File
import kotlin.math.abs
import kotlin.math.max
import kotlin.test.Test
import kotlin.test.assertTrue

/**
 * Диагностика обусловленности матриц аппроксимационного соотношения M_k и разбиения
 * единицы при измельчении равномерной сетки. Результат — TSV `build/reports/conditioning-diagnostic.tsv`
 * (колонки `sys n k cond pu`), сводка печатается одной строкой на конфигурацию.
 */
@Tag("fast")
class ConditioningDiagnosticTest {

    @Test
    fun conditioningAndPartitionOfUnityOnRefinement() {
        val backend = Backends.default()
        val lines = ArrayList<String>()
        lines += "sys\tn\tk\tcond\tpu"
        val buildFailures = ArrayList<String>()
        val puFailures = ArrayList<String>()
        val configs = ArrayList<Pair<GeneratingSystem, Grid>>()
        for (sys in listOf(GeneratingSystem.B, GeneratingSystem.H)) {
            for (n in listOf(10, 100, 1000, 10000)) configs += sys to Grid.uniform(n, 0.0, 1.0)
        }
        configs += GeneratingSystem.B to Grid.uniform(100, 100.0, 101.0)
        for ((sys, grid) in configs) {
            val n = grid.n
            val a = grid.x(0); val b = grid.x(n)
            val label = if (a == 0.0) sys.name else "${sys.name}[$a,$b]"
            val basis = try { MinimalSplineBasis(sys, grid) } catch (e: Exception) {
                val msg = "${e::class.simpleName}: ${e.message}"
                lines += "$label\t$n\t-\t-\t-\t$msg"
                buildFailures += "$label n=$n: $msg"
                println("conditioning-diagnostic: $label n=$n НЕ ПОСТРОЕН: $msg")
                continue
            }
            val ones = DoubleArray(n + 2) { 1.0 }
            var pu = 0.0
            for (i in 1..50) {
                val t = a + (b - a) * i / 51.0
                pu = max(pu, abs(basis.evalSpline(ones, t) - 1.0))
            }
            val conds = ArrayList<Double>()
            for (k in listOf(0, n / 2, n - 1)) {
                val cond = Conditioning.conditionEstimate(basis.approximationMatrix(k), backend = backend).condInf
                conds += cond
                lines += "$label\t$n\t$k\t${"%.3e".format(cond)}\t${"%.3e".format(pu)}"
            }
            println("conditioning-diagnostic: $label n=$n cond(k=0,n/2,n-1)=${conds.map { "%.2e".format(it) }} pu=${"%.2e".format(pu)}")
            if (pu > 1e-6) puFailures += "$label n=$n: pu=$pu"
        }
        val out = File("build/reports/conditioning-diagnostic.tsv")
        out.parentFile.mkdirs()
        out.writeText(lines.joinToString("\n") + "\n")
        assertTrue(puFailures.isEmpty(), "разбиение единицы нарушено: $puFailures")
        if (buildFailures.isNotEmpty()) println("conditioning-diagnostic: не построены: $buildFailures")
        // 1c: критерий достоверности — оценка числа обусловленности M_k, не превышающая
        // MinimalSplineBasis.MAX_CONDITION. B на [100,101] отвергается с исключением о числе
        // обусловленности (глобальные координаты порождающей системы); на [0,1] базис строится.
        // Измерено (1c): cond(M_k) на [0,1] растёт как n²; при n = 10⁴ она равна ~3·10⁸ (B) и ~10⁹ (H),
        // то есть превышает MAX_CONDITION = 10⁸, и базис в глобальных координатах не строится.
        // Ожидается построение после перехода к локальным координатам (этап 1d); здесь фиксируется
        // только то, что любой отказ вызван именно критерием обусловленности.
        assertTrue(
            buildFailures.all { it.contains("число обусловленности") },
            "отказ построения не по числу обусловленности: $buildFailures",
        )
        assertTrue(
            buildFailures.any { it.startsWith("B[100.0,101.0] n=100") },
            "B на [100,101] должен отвергаться по числу обусловленности: $buildFailures",
        )
        val smallN = buildFailures.filter { f -> listOf(10, 100, 1000).any { f.contains(" n=$it:") } && !f.startsWith("B[") }
        assertTrue(smallN.isEmpty(), "базис на [0,1] при n ≤ 10³ должен строиться: $smallN")
    }
}
