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

/**
 * algebra.json: nonDegenerate — bit for bit (integer logic); phi/wronskian — relative
 * tolerance [GoldenCompare.TRANSCENDENTAL_TOLERANCE] (the T system is computed via `Math.sin`/`Math.cos`).
 */
@Tag("fast")
class GoldenAlgebraTest {
    @TestFactory
    fun algebra(): List<DynamicTest> {
        val root = readGolden("algebra.json")
        val expNd = section(root, "nonDegenerate")
        val gotNd = GoldenCompute.nonDegenerateCases()
        val nd = expNd.keys.map { key -> dynamicTest("nonDegenerate/$key") { compare(expNd[key], gotNd[key], key, Mode.BITS) } }
        val expPhi = section(root, "phi")
        val gotPhi = GoldenCompute.phiCases()
        val phi = expPhi.keys.map { key -> dynamicTest("phi/$key") { compare(expPhi[key], gotPhi[key], key, Mode.REL, TRANSCENDENTAL_TOLERANCE) } }
        return nd + phi + dynamicTest("case set") {
            compare(expNd.keys.sorted(), gotNd.keys.sorted(), "nonDegenerate", Mode.BITS)
            compare(expPhi.keys.sorted(), gotPhi.keys.sorted(), "phi", Mode.BITS)
        }
    }
}
