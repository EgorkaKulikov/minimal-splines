package splines.golden

import org.junit.jupiter.api.DynamicTest
import org.junit.jupiter.api.DynamicTest.dynamicTest
import org.junit.jupiter.api.Tag
import org.junit.jupiter.api.TestFactory
import splines.golden.GoldenCompare.Mode
import splines.golden.GoldenCompare.compare
import splines.golden.GoldenIo.readGolden
import splines.golden.GoldenIo.section

/** metrics.json: errorEh — относительно 1e-12; orders и constCh — побитово (явные формулы). */
@Tag("fast")
class GoldenMetricsTest {
    @TestFactory
    fun metrics(): List<DynamicTest> {
        val root = readGolden("metrics.json")
        val expEh = section(root, "errorEh")
        val got by lazy { GoldenCompute.errorEhCases() }
        val eh = expEh.keys.map { key -> dynamicTest("errorEh/$key") { compare(expEh[key], got[key], key, Mode.REL) } }
        return eh + listOf(
            dynamicTest("errorEh/набор случаев") { compare(expEh.keys.sorted(), got.keys.sorted(), "errorEh", Mode.BITS) },
            dynamicTest("orders") { compare(section(root, "orders"), GoldenCompute.ordersCases(), "orders", Mode.BITS) },
            dynamicTest("constCh") { compare(section(root, "constCh"), GoldenCompute.constChCases(), "constCh", Mode.BITS) },
        )
    }
}
