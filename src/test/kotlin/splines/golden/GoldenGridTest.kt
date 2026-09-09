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
import kotlin.test.fail

/** grid.json: узлы, шаг, признаки совпадения — побитово (целочисленная логика и явные формулы). */
@Tag("fast")
class GoldenGridTest {
    @TestFactory
    fun grids(): List<DynamicTest> {
        val exp = section(readGolden("grid.json"), "grids")
        val cases = exp.keys.map { key ->
            dynamicTest(key) {
                val g = GoldenInputs.grids[key] ?: fail("сетка $key отсутствует во входах")
                compare(exp[key], GoldenCompute.gridCase(g), key, Mode.BITS)
            }
        }
        return cases + dynamicTest("набор сеток") { assertEquals(GoldenInputs.grids.keys, exp.keys) }
    }
}
