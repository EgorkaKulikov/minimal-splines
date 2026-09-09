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
 * metrics.json: errorEh — относительно [GoldenCompare.BASIS_TOLERANCE] (величина проходит через значения
 * сплайна и наследует погрешность эталона базиса); orders и constCh — побитово (явные формулы).
 */
@Tag("fast")
class GoldenMetricsTest {
    @TestFactory
    fun metrics(): List<DynamicTest> {
        val root = readGolden("metrics.json")
        val expEh = section(root, "errorEh")
        val got by lazy { GoldenCompute.errorEhCases() }
        // errorEh вычисляется через значения сплайна и наследует погрешность эталона базиса
        // (формулы Крамера в 0.1.0), поэтому сравнивается с тем же допуском, что и basis.
        val eh = expEh.keys.map { key ->
            dynamicTest("errorEh/$key") { compare(expEh[key], got[key], key, Mode.REL, GoldenCompare.BASIS_TOLERANCE) }
        }
        return eh + listOf(
            dynamicTest("errorEh/набор случаев") { compare(expEh.keys.sorted(), got.keys.sorted(), "errorEh", Mode.BITS) },
            dynamicTest("orders") { compare(section(root, "orders"), GoldenCompute.ordersCases(), "orders", Mode.BITS) },
            dynamicTest("constCh") { compare(section(root, "constCh"), GoldenCompute.constChCases(), "constCh", Mode.BITS) },
        )
    }
}
