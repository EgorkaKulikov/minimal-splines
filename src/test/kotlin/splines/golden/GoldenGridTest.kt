package splines.golden

import org.junit.jupiter.api.DynamicTest
import org.junit.jupiter.api.DynamicTest.dynamicTest
import org.junit.jupiter.api.Tag
import org.junit.jupiter.api.TestFactory
import splines.golden.GoldenCompare.Mode
import splines.golden.GoldenCompare.TRANSCENDENTAL_TOLERANCE
import splines.golden.GoldenCompare.compare
import splines.golden.GoldenIo.readGolden
import splines.golden.GoldenIo.section
import kotlin.test.assertEquals
import kotlin.test.fail

/**
 * grid.json: nodes, step, coincidence flags — bit for bit (integer logic and explicit formulas);
 * the `quasiUniform` and `geometric` grids (nodes computed via `Math.sin`/`Math.pow`) — relative tolerance
 * [GoldenCompare.TRANSCENDENTAL_TOLERANCE].
 */
@Tag("fast")
class GoldenGridTest {
    @TestFactory
    fun grids(): List<DynamicTest> {
        val exp = section(readGolden("grid.json"), "grids")
        val cases = exp.keys.map { key ->
            dynamicTest(key) {
                val g = GoldenInputs.grids[key] ?: fail("grid $key is missing from the inputs")
                if (key.startsWith("quasiUniform-") || key.startsWith("geometric-")) {
                    compare(exp[key], GoldenCompute.gridCase(g), key, Mode.REL, TRANSCENDENTAL_TOLERANCE)
                } else {
                    compare(exp[key], GoldenCompute.gridCase(g), key, Mode.BITS)
                }
            }
        }
        return cases + dynamicTest("grid set") { assertEquals(GoldenInputs.grids.keys, exp.keys) }
    }
}
