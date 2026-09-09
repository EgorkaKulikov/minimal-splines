package splines.golden

import org.junit.jupiter.api.DynamicTest
import org.junit.jupiter.api.DynamicTest.dynamicTest
import org.junit.jupiter.api.Tag
import org.junit.jupiter.api.TestFactory
import splines.golden.GoldenCompare.Mode
import splines.golden.GoldenCompare.compare
import splines.golden.GoldenIo.readGolden
import splines.golden.GoldenIo.section
import kotlin.test.assertEquals

/**
 * basis.json: значения базиса и сплайна — относительно [GoldenCompare.BASIS_TOLERANCE]
 * (коэффициенты проходят через обращение 3x3 матриц; эталон 0.1.0 получен по формулам Крамера,
 * см. обоснование допуска); `interval` — точно.
 */
@Tag("fast")
class GoldenBasisTest {
    @TestFactory
    fun basis(): List<DynamicTest> {
        val exp = section(readGolden("basis.json"), "cases")
        val inputs = GoldenInputs.basisCases().associate { (s, gk, g) -> "$s-$gk" to (GoldenInputs.systems.getValue(s) to g) }
        val cases = exp.keys.map { key ->
            dynamicTest(key) {
                val (sys, g) = inputs.getValue(key)
                compare(exp[key], GoldenCompute.basisCase(sys, g), key, Mode.REL, GoldenCompare.BASIS_TOLERANCE)
            }
        }
        return cases + dynamicTest("набор случаев") { assertEquals(inputs.keys, exp.keys) }
    }
}
