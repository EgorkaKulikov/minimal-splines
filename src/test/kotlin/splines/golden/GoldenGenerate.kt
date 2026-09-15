package splines.golden

import org.junit.jupiter.api.Tag
import org.junit.jupiter.api.Test
import splines.golden.GoldenCompute.basisCase
import splines.golden.GoldenCompute.functionalsCase
import splines.golden.GoldenCompute.gridCase
import splines.golden.GoldenIo.writeGolden
import java.io.File

/**
 * Generation of the golden references: `./gradlew regenerateGolden`. Not part of `test`
 * (the `golden-generate` tag is excluded); the directory is set by the `golden.dir` property.
 */
@Tag("golden-generate")
class GoldenGenerate {
    private val dir = File(System.getProperty("golden.dir") ?: error("system property golden.dir is not set"))

    private fun root(vararg sections: Pair<String, Any?>): Map<String, Any?> =
        linkedMapOf<String, Any?>("generatedWith" to GoldenCompute.meta()).apply { putAll(sections) }

    @Test
    fun grid() = writeGolden(dir, "grid.json", root(
        "grids" to GoldenInputs.grids.mapValues { gridCase(it.value) },
    ))

    @Test
    fun basis() = writeGolden(dir, "basis.json", root(
        "cases" to GoldenInputs.basisCases().associate { (s, gk, g) -> "$s-$gk" to basisCase(GoldenInputs.systems.getValue(s), g) },
    ))

    @Test
    fun functionals() = writeGolden(dir, "functionals.json", root(
        "cases" to GoldenInputs.basisCases(16).associate { (s, gk, g) -> "$s-$gk" to functionalsCase(GoldenInputs.systems.getValue(s), g) },
    ))

    @Test
    fun metrics() = writeGolden(dir, "metrics.json", root(
        "errorEh" to GoldenCompute.errorEhCases(),
    ))

    @Test
    fun algebra() = writeGolden(dir, "algebra.json", root(
        "nonDegenerate" to GoldenCompute.nonDegenerateCases(),
        "phi" to GoldenCompute.phiCases(),
    ))
}
