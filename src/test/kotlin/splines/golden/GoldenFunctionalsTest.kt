package splines.golden

import org.junit.jupiter.api.DynamicTest
import org.junit.jupiter.api.DynamicTest.dynamicTest
import org.junit.jupiter.api.Tag
import org.junit.jupiter.api.TestFactory
import splines.MinimalSplineBasis
import splines.golden.GoldenCompare.Mode
import splines.golden.GoldenCompare.compare
import splines.golden.GoldenIo.obj
import splines.golden.GoldenIo.readGolden
import splines.golden.GoldenIo.section
import kotlin.test.assertEquals

/**
 * functionals.json: коэффициенты проекторов, cChi, closedFormInternal — относительно
 * [GoldenCompare.BASIS_TOLERANCE] (семейства решают малые СЛАУ бэкендом numerical-core, эталон 0.1.0
 * получен по формулам Крамера, см. обоснование допуска); флаги и классы исключений — точно.
 */
@Tag("fast")
class GoldenFunctionalsTest {
    @TestFactory
    fun functionals(): List<DynamicTest> {
        val exp = section(readGolden("functionals.json"), "cases")
        val inputs = GoldenInputs.basisCases(16).associate { (s, gk, g) -> "$s-$gk" to (GoldenInputs.systems.getValue(s) to g) }
        val cases = exp.keys.flatMap { key ->
            val fams = obj(exp[key])
            fams.keys.map { fam ->
                dynamicTest("$key/$fam") {
                    val (sys, g) = inputs.getValue(key)
                    compare(
                        fams[fam], GoldenCompute.familyCase(fam, MinimalSplineBasis(sys, g)), "$key/$fam", Mode.REL,
                        GoldenCompare.BASIS_TOLERANCE,
                    )
                }
            }
        }
        return cases + dynamicTest("набор случаев") {
            assertEquals(inputs.keys, exp.keys)
            for (key in exp.keys) assertEquals(GoldenCompute.familyKeys.toSet(), obj(exp[key]).keys, key)
        }
    }
}
