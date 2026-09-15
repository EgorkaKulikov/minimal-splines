package splines.golden

import org.junit.jupiter.api.DynamicTest
import org.junit.jupiter.api.DynamicTest.dynamicTest
import org.junit.jupiter.api.Tag
import org.junit.jupiter.api.TestFactory
import splines.golden.GoldenCompare.Mode
import splines.golden.GoldenCompare.compare
import splines.golden.GoldenIo.readGolden
import splines.golden.GoldenIo.section

/**
 * metrics.json: errorEh — compared against [GoldenCompare.BASIS_TOLERANCE] (the quantity goes through the
 * spline values and inherits the error of the basis golden reference).
 */
@Tag("fast")
class GoldenMetricsTest {
    @TestFactory
    fun metrics(): List<DynamicTest> {
        val root = readGolden("metrics.json")
        val expEh = section(root, "errorEh")
        val got by lazy { GoldenCompute.errorEhCases() }
        // errorEh is computed through the spline values and inherits the error of the basis golden reference
        // (Cramer's rule in 0.1.0), so it is compared with the same tolerance as basis.
        val eh = expEh.keys.map { key ->
            dynamicTest("errorEh/$key") { compare(expEh[key], got[key], key, Mode.REL, GoldenCompare.BASIS_TOLERANCE) }
        }
        return eh + listOf(
            dynamicTest("errorEh/case set") { compare(expEh.keys.sorted(), got.keys.sorted(), "errorEh", Mode.BITS) },
        )
    }
}
