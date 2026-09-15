package splines

import numerics.NumericsContext
import numerics.backend.Backends
import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

/**
 * The default BLAS/LAPACK implementation matches the one requested via the `numerics.backend`
 * property. In CI, for the native variant this test runs as a separate step before the main test
 * run: without the system libraries netlib falls back to F2J, and the native variant would then
 * check exactly the same thing as the java one.
 */
@Tag("fast")
class BackendSelectionTest {
    @Test
    fun defaultMatchesRequestedProperty() {
        val requested = System.getProperty("numerics.backend")?.trim()?.lowercase()
        val actual = Backends.default()
        val description = Backends.describe()
        when (requested) {
            "java" -> assertFalse(actual.isNative, description)
            "native" -> assertTrue(actual.isNative, description)
            else -> assertEquals(Backends.isNativeAvailable(), actual.isNative, description)
        }
        assertEquals(actual, NumericsContext.default().backend)
        assertTrue(description.contains(actual.name), description)
    }
}
