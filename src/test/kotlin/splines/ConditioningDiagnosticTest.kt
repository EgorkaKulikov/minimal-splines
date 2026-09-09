package splines

import numerics.Conditioning
import numerics.DenseMatrix
import numerics.backend.Backends
import org.junit.jupiter.api.Tag
import splines.functionals.ProjFunctionals
import java.io.File
import kotlin.math.abs
import kotlin.math.cosh
import kotlin.math.max
import kotlin.math.sinh
import kotlin.test.Test
import kotlin.test.assertTrue

/**
 * Обусловленность матриц аппроксимационного соотношения в локальных координатах интервала
 * (T_k M_k) при измельчении сетки и при изменении положения и масштаба отрезка. Для сравнения
 * записывается и число обусловленности глобальной матрицы M_k. Результат — TSV
 * `build/reports/conditioning-diagnostic.tsv` (колонки `sys a b n condLocalMax condGlobal(k=0,n/2,n-1) pu`).
 */
@Tag("fast")
class ConditioningDiagnosticTest {
    private val backend = Backends.default()

    private fun cond(m: DenseMatrix): Double = Conditioning.conditionEstimate(m, backend = backend).condInf

    private fun partitionOfUnityDefect(basis: MinimalSplineBasis): Double {
        val grid = basis.grid
        val ones = DoubleArray(grid.n + 2) { 1.0 }
        var pu = 0.0
        for (i in 1..50) {
            val t = grid.a + (grid.b - grid.a) * i / 51.0
            pu = max(pu, abs(basis.evalSpline(ones, t) - 1.0))
        }
        return pu
    }

    private fun maxLocalCond(basis: MinimalSplineBasis): Double =
        (0 until basis.n).maxOf { cond(basis.localApproximationMatrix(it)) }

    private fun globalConds(basis: MinimalSplineBasis): String =
        listOf(0, basis.n / 2, basis.n - 1).joinToString(",") { k ->
            try { "%.2e".format(cond(basis.approximationMatrix(k))) } catch (e: Exception) { "-" }
        }

    private fun record(lines: MutableList<String>, sys: GeneratingSystem, grid: Grid): Pair<Double, Double> {
        val basis = MinimalSplineBasis(sys, grid)
        val condLocal = maxLocalCond(basis)
        val pu = partitionOfUnityDefect(basis)
        lines += "${sys.name}\t${grid.a}\t${grid.b}\t${grid.n}\t${"%.3e".format(condLocal)}\t${globalConds(basis)}\t${"%.3e".format(pu)}"
        return condLocal to pu
    }

    private fun write(name: String, lines: List<String>) {
        val out = File("build/reports/$name")
        out.parentFile.mkdirs()
        out.writeText(lines.joinToString("\n") + "\n")
    }

    /** cond(T_k M_k) на [0,1] не превышает 10³ и не растёт при измельчении: cond(n=10⁴)/cond(n=10) ≤ 10. */
    @Test
    fun localConditioningIsIndependentOfRefinement() {
        val lines = arrayListOf("sys\ta\tb\tn\tcondLocalMax\tcondGlobal\tpu")
        val failures = ArrayList<String>()
        for (sys in listOf(GeneratingSystem.B, GeneratingSystem.H, GeneratingSystem.T)) {
            val conds = LinkedHashMap<Int, Double>()
            for (n in listOf(10, 100, 1000, 10000)) {
                val (condLocal, pu) = record(lines, sys, Grid.uniform(n, 0.0, 1.0))
                conds[n] = condLocal
                if (condLocal > 1e3) failures += "${sys.name} n=$n: cond=$condLocal"
                if (pu > 1e-9) failures += "${sys.name} n=$n: разбиение единицы pu=$pu"
            }
            val ratio = conds.getValue(10000) / conds.getValue(10)
            if (ratio > 10.0) failures += "${sys.name}: cond растёт с n, cond(10⁴)/cond(10)=$ratio"
        }
        write("conditioning-diagnostic.tsv", lines)
        assertTrue(failures.isEmpty(), "обусловленность в локальных координатах: $failures")
    }

    /** cond(T_k M_k) не зависит от положения и масштаба отрезка: B на [100,101], [0,10⁶], [0,10⁻⁶]; H на [100,101]. */
    @Test
    fun localConditioningIsIndependentOfSegment() {
        val lines = arrayListOf("sys\ta\tb\tn\tcondLocalMax\tcondGlobal\tpu")
        val failures = ArrayList<String>()
        for ((a, b) in listOf(100.0 to 101.0, 0.0 to 1e6, 0.0 to 1e-6)) {
            val (condLocal, pu) = record(lines, GeneratingSystem.B, Grid.uniform(100, a, b))
            if (condLocal > 1e3) failures += "B[$a,$b]: cond=$condLocal"
            if (pu > 1e-9) failures += "B[$a,$b]: разбиение единицы pu=$pu"
        }
        val (condH, puH) = record(lines, GeneratingSystem.H, Grid.uniform(100, 100.0, 101.0))
        if (condH > 1e3) failures += "H[100,101]: cond=$condH"
        if (puH > 1e-8) failures += "H[100,101]: разбиение единицы pu=$puH"
        write("conditioning-diagnostic-segments.tsv", lines)
        assertTrue(failures.isEmpty(), "обусловленность в локальных координатах: $failures")
    }

    /** При n = 10⁴ разбиение единицы и воспроизведение элемента span phi выполняются с точностью 10⁻⁹. */
    @Test
    fun fineGridReproducesGeneratingSpan() {
        val cases = listOf(
            Triple(GeneratingSystem.B, { t: Double -> t * t }, { t: Double -> 2.0 * t } to { _: Double -> 2.0 }),
            Triple(GeneratingSystem.H, { t: Double -> cosh(t) }, { t: Double -> sinh(t) } to { t: Double -> cosh(t) }),
        )
        for ((sys, f, derivs) in cases) {
            val grid = Grid.uniform(10000, 0.0, 1.0)
            val basis = MinimalSplineBasis(sys, grid)
            val pu = partitionOfUnityDefect(basis)
            assertTrue(pu <= 1e-9, "${sys.name} n=10⁴: разбиение единицы pu=$pu")
            val c = ProjFunctionals(basis).projectorCoeffs(f, derivs.first, derivs.second)
            var err = 0.0
            for (i in 0..50) {
                val t = i / 50.0
                err = max(err, abs(basis.evalSpline(c, t) - f(t)))
            }
            assertTrue(err <= 1e-9, "${sys.name} n=10⁴: погрешность на span phi $err")
        }
    }
}
