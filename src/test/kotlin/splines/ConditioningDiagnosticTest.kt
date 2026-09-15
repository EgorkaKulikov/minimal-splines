package splines

import numerics.Conditioning
import numerics.DenseMatrix
import numerics.backend.Backends
import org.junit.jupiter.api.Tag
import splines.functionals.AveragingFunctionals
import splines.functionals.DeBoorFixFunctionals
import splines.functionals.DiscreteDeBoorFixFunctionals
import splines.functionals.FunctionalFamily
import splines.functionals.ProjFunctionals
import splines.functionals.ThreePointFunctionals
import splines.metrics.errorEh
import java.io.File
import kotlin.math.abs
import kotlin.math.cosh
import kotlin.math.max
import kotlin.math.sinh
import kotlin.test.Test
import kotlin.test.assertTrue

/**
 * Conditioning of the approximation relation matrices in the local coordinates of the interval
 * (T_k M_k) under grid refinement and under changes of the position and scale of the interval. For
 * comparison, the condition number of the global matrix M_k is recorded as well. The result is the
 * TSV `build/reports/conditioning-diagnostic.tsv` (columns `sys a b n condLocalMax condGlobal(k=0,n/2,n-1) pu`).
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

    /** cond(T_k M_k) on [0,1] does not exceed 10³ and does not grow under refinement: cond(n=10⁴)/cond(n=10) ≤ 10. */
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
                if (pu > 1e-9) failures += "${sys.name} n=$n: partition of unity pu=$pu"
            }
            val ratio = conds.getValue(10000) / conds.getValue(10)
            if (ratio > 10.0) failures += "${sys.name}: cond grows with n, cond(10⁴)/cond(10)=$ratio"
        }
        write("conditioning-diagnostic.tsv", lines)
        assertTrue(failures.isEmpty(), "conditioning in local coordinates: $failures")
    }

    /** cond(T_k M_k) does not depend on the position and scale of the interval: B on [100,101], [0,10⁶], [0,10⁻⁶]; H on [100,101]. */
    @Test
    fun localConditioningIsIndependentOfSegment() {
        val lines = arrayListOf("sys\ta\tb\tn\tcondLocalMax\tcondGlobal\tpu")
        val failures = ArrayList<String>()
        for ((a, b) in listOf(100.0 to 101.0, 0.0 to 1e6, 0.0 to 1e-6)) {
            val (condLocal, pu) = record(lines, GeneratingSystem.B, Grid.uniform(100, a, b))
            if (condLocal > 1e3) failures += "B[$a,$b]: cond=$condLocal"
            if (pu > 1e-9) failures += "B[$a,$b]: partition of unity pu=$pu"
        }
        val (condH, puH) = record(lines, GeneratingSystem.H, Grid.uniform(100, 100.0, 101.0))
        if (condH > 1e3) failures += "H[100,101]: cond=$condH"
        if (puH > 1e-8) failures += "H[100,101]: partition of unity pu=$puH"
        write("conditioning-diagnostic-segments.tsv", lines)
        assertTrue(failures.isEmpty(), "conditioning in local coordinates: $failures")
    }

    /** For n = 10⁴ the partition of unity and the reproduction of an element of span phi hold to 10⁻⁹. */
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
            assertTrue(pu <= 1e-9, "${sys.name} n=10⁴: partition of unity pu=$pu")
            val c = ProjFunctionals(basis).projectorCoeffs(f, derivs.first, derivs.second)
            var err = 0.0
            for (i in 0..50) {
                val t = i / 50.0
                err = max(err, abs(basis.evalSpline(c, t) - f(t)))
            }
            assertTrue(err <= 1e-9, "${sys.name} n=10⁴: error on span phi $err")
        }
    }

    /**
     * All five families of functionals are built for B on [0,1] with n = 10⁴ and for B, H, T on
     * [100,101] with n = 100,
     * and every quasi-projector reproduces f = 1 + 2 rho + 3 sigma to within 10⁻⁹ max|f|
     * in the errorEh metric. For xitilde the accuracy on span phi is checked only for B and outside
     * the boundary layer of three steps: on a uniform grid the central difference reproduces the
     * derivative of a quadratic polynomial exactly, while for H and T it gives an error of order h³ by
     * the construction of the family (the order is checked in ConvergenceOrderTest); near the edges the
     * derivative is replaced by a one-sided difference over multiple nodes, and the error on span phi
     * there is of order h². For H and T only the construction of xitilde is checked, and the error is
     * written to the report.
     * The result is the TSV `build/reports/functionals-diagnostic.tsv`.
     */
    @Test
    fun functionalsBuildOnFineGridsAndShiftedIntervals() {
        val lines = arrayListOf("system\ta\tb\tn\tfamily\tstatus\trelErr\trelErrInterior")
        val failures = ArrayList<String>()
        val cases = listOf(GeneratingSystem.B to Grid.uniform(10000, 0.0, 1.0)) +
            listOf(GeneratingSystem.B, GeneratingSystem.H, GeneratingSystem.T).map { it to Grid.uniform(100, 100.0, 101.0) }
        for ((sys, grid) in cases) {
            val basis = MinimalSplineBasis(sys, grid)
            val f = { t: Double -> sys.phi(t).let { it[0] + 2.0 * it[1] + 3.0 * it[2] } }
            val fD = { t: Double -> sys.phiD(t).let { 2.0 * it[1] + 3.0 * it[2] } }
            val fDD = { t: Double -> sys.phiDD(t).let { 2.0 * it[1] + 3.0 * it[2] } }
            var fMax = 0.0
            for (i in 0..100) fMax = max(fMax, abs(f(grid.a + (grid.b - grid.a) * i / 100.0)))
            val families: List<Pair<String, () -> FunctionalFamily>> = listOf(
                "theta" to { ProjFunctionals(basis) },
                "xi" to { DeBoorFixFunctionals(basis, r = 1) },
                "xitilde" to { DiscreteDeBoorFixFunctionals(basis, r = 1) },
                "mu" to { AveragingFunctionals(basis) },
                "lambda" to { ThreePointFunctionals(basis) },
            )
            for ((name, make) in families) {
                val tag = "${sys.name}[${grid.a},${grid.b}] n=${grid.n} $name"
                val family = try { make() } catch (e: Exception) {
                    lines += "${sys.name}\t${grid.a}\t${grid.b}\t${grid.n}\t$name\tfailed: ${e.message}\t-\t-"
                    failures += "$tag: is not built (${e.message})"
                    continue
                }
                val c = family.projectorCoeffs(f, fD, fDD)
                val relErr = errorEh(f, { t -> basis.evalSpline(c, t) }, grid, refinement = 10) / fMax
                val h = (grid.b - grid.a) / grid.n
                var interior = 0.0
                for (i in 0..(10 * grid.n)) {
                    val t = grid.a + (grid.b - grid.a) * i / (10 * grid.n)
                    if (t >= grid.a + 3 * h && t <= grid.b - 3 * h) interior = max(interior, abs(basis.evalSpline(c, t) - f(t)))
                }
                val relErrInterior = interior / fMax
                lines += "${sys.name}\t${grid.a}\t${grid.b}\t${grid.n}\t$name\tok\t${"%.3e".format(relErr)}\t${"%.3e".format(relErrInterior)}"
                val checked = if (name == "xitilde") relErrInterior else relErr
                val exactOnSpan = name != "xitilde" || sys === GeneratingSystem.B
                if (exactOnSpan && checked > 1e-9) failures += "$tag: relative error on span phi $checked"
            }
        }
        write("functionals-diagnostic.tsv", lines)
        assertTrue(failures.isEmpty(), "families of functionals on a fine grid and a shifted interval: $failures")
    }
}
