package splines

import org.junit.jupiter.api.Tag
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertTrue

@Tag("fast")
class GeneratingSystemTest {
    @Test fun wronskianNonZero() {
        for (sys in listOf(GeneratingSystem.B, GeneratingSystem.H, GeneratingSystem.T)) {
            for (t in listOf(0.1, 0.3, 0.5, 0.7, 0.9)) {
                assertTrue(abs(sys.wronskian(t)) > 1e-9, "wronskian zero for ${sys.name} at $t")
            }
        }
    }
}
