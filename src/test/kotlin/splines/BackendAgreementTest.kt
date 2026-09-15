package splines

import numerics.NumericsContext
import numerics.backend.Backends
import org.junit.jupiter.api.Assumptions.assumeTrue
import org.junit.jupiter.api.DynamicTest
import org.junit.jupiter.api.DynamicTest.dynamicTest
import org.junit.jupiter.api.Tag
import org.junit.jupiter.api.TestFactory
import splines.functionals.AveragingFunctionals
import splines.functionals.DeBoorFixFunctionals
import splines.functionals.DiscreteDeBoorFixFunctionals
import splines.functionals.FunctionalFamily
import splines.functionals.ProjFunctionals
import splines.functionals.ThreePointFunctionals
import java.io.File
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.exp
import kotlin.math.max
import kotlin.math.sin
import kotlin.test.assertTrue

/**
 * Agreement between the native (Accelerate/OpenBLAS) and the pure Java BLAS/LAPACK
 * implementations: the basis ω_j and the coefficients of the five functional families, built with
 * the two contexts, must coincide up to rounding errors. The tests are skipped if the native
 * library is unavailable. The measured maximal discrepancies are written to
 * `build/reports/backend-agreement.tsv` (columns `sys grid kind name maxRelDiff`).
 */
@Tag("fast")
class BackendAgreementTest {
    private val systems = listOf(GeneratingSystem.B, GeneratingSystem.H, GeneratingSystem.T)
    private val grids = listOf("uniform" to Grid.uniform(32, 0.0, 1.0), "graded" to Grid.graded(32, 0.0, 1.0))

    private val f: (Double) -> Double = { t -> exp(sin(3.0 * t)) }
    private val fD: (Double) -> Double = { t -> 3.0 * cos(3.0 * t) * exp(sin(3.0 * t)) }
    private val fDD: (Double) -> Double = { t ->
        val c = cos(3.0 * t)
        (9.0 * c * c - 9.0 * sin(3.0 * t)) * exp(sin(3.0 * t))
    }

    private val basisTol = 1e-12
    private val coeffTol = 1e-10

    private companion object {
        /** The file is recreated once per JVM run and afterwards only appended to. */
        private val reportFile: File by lazy {
            File("build/reports/backend-agreement.tsv").also { file ->
                file.parentFile.mkdirs()
                file.writeText("sys\tgrid\tkind\tname\tmaxRelDiff\n")
            }
        }
    }

    private fun report(): File = reportFile

    private fun relDiff(a: DoubleArray, b: DoubleArray): Double {
        require(a.size == b.size)
        var scale = 0.0
        var diff = 0.0
        for (i in a.indices) {
            scale = max(scale, abs(b[i]))
            diff = max(diff, abs(a[i] - b[i]))
        }
        return if (scale == 0.0) diff else diff / scale
    }

    private fun families(basis: MinimalSplineBasis, ctx: NumericsContext): List<FunctionalFamily> = listOf(
        ProjFunctionals(basis, ctx),
        DeBoorFixFunctionals(basis, 1, ctx),
        DiscreteDeBoorFixFunctionals(basis, 1, ctx),
        AveragingFunctionals(basis, 0.5, ctx),
        ThreePointFunctionals(basis, 0.5, ctx),
    )

    @TestFactory
    fun nativeAndJavaAgree(): List<DynamicTest> = systems.flatMap { sys ->
        grids.map { (gridName, grid) ->
            dynamicTest("${sys.name} / $gridName: basis and five families") {
                assumeTrue(Backends.isNativeAvailable(), "native BLAS/LAPACK implementation is unavailable")
                val ctxN = NumericsContext(backend = Backends.native())
                val ctxJ = NumericsContext(backend = Backends.java())
                val basisN = MinimalSplineBasis(sys, grid, ctxN)
                val basisJ = MinimalSplineBasis(sys, grid, ctxJ)
                val out = report()

                // Basis: ω_j at 100 points of the interval, for all j.
                val points = DoubleArray(100) { grid.a + (grid.b - grid.a) * (it + 0.5) / 100.0 }
                var basisDiff = 0.0
                for (j in -2..grid.n - 1) {
                    val wN = DoubleArray(points.size) { basisN.omega(j, points[it]) }
                    val wJ = DoubleArray(points.size) { basisJ.omega(j, points[it]) }
                    basisDiff = max(basisDiff, relDiff(wN, wJ))
                }
                out.appendText("${sys.name}\t$gridName\tbasis\tomega\t$basisDiff\n")
                assertTrue(basisDiff <= basisTol, "${sys.name}/$gridName: basis discrepancy $basisDiff > $basisTol")

                // Functionals: quasi-projector coefficients for f = exp(sin 3t).
                val famN = families(basisN, ctxN)
                val famJ = families(basisJ, ctxJ)
                for ((fN, fJ) in famN.zip(famJ)) {
                    val cN = fN.projectorCoeffs(f, fD, fDD)
                    val cJ = fJ.projectorCoeffs(f, fD, fDD)
                    val d = relDiff(cN, cJ)
                    out.appendText("${sys.name}\t$gridName\tfunctional\t${fN.name}\t$d\n")
                    assertTrue(d <= coeffTol, "${sys.name}/$gridName/${fN.name}: coefficient discrepancy $d > $coeffTol")
                }
            }
        }
    }
}
