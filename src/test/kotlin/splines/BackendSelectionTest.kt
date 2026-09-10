package splines

import numerics.NumericsContext
import numerics.backend.Backends
import org.junit.jupiter.api.Tag
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

/**
 * Реализация BLAS/LAPACK по умолчанию соответствует запрошенной свойством `numerics.backend`.
 * В CI на варианте native этот тест выполняется отдельным шагом до основного запуска тестов: без
 * системных библиотек netlib переключается на F2J, и вариант native проверял бы то же, что и java.
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
