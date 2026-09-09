package splines

// ============================================================================
// 3. ПОРОЖДАЮЩАЯ ВЕКТОР-ФУНКЦИЯ phi(t) = (1, rho(t), sigma(t))^T
// ============================================================================

/**
 * Порождающая вектор-функция phi(t) = (1, rho(t), sigma(t))^T и её производные
 * до второго порядка. phi_0 == 1 обеспечивает разбиение единицы.
 */
class GeneratingSystem(
    val name: String,
    val rho: (Double) -> Double,
    val sigma: (Double) -> Double,
    val rhoD: (Double) -> Double,
    val sigmaD: (Double) -> Double,
    val rhoDD: (Double) -> Double,
    val sigmaDD: (Double) -> Double,
) {
    /** phi(t) = (1, rho(t), sigma(t)). */
    fun phi(t: Double): DoubleArray = doubleArrayOf(1.0, rho(t), sigma(t))

    /** phi'(t) = (0, rho'(t), sigma'(t)). */
    fun phiD(t: Double): DoubleArray = doubleArrayOf(0.0, rhoD(t), sigmaD(t))

    /** phi''(t) = (0, rho''(t), sigma''(t)). */
    fun phiDD(t: Double): DoubleArray = doubleArrayOf(0.0, rhoDD(t), sigmaDD(t))

    /** Вронскиан det(phi, phi', phi'') — проверка невырожденности. */
    fun wronskian(t: Double): Double {
        val a = phi(t); val b = phiD(t); val c = phiDD(t)
        return a[0] * (b[1] * c[2] - b[2] * c[1]) +
            a[1] * (b[2] * c[0] - b[0] * c[2]) +
            a[2] * (b[0] * c[1] - b[1] * c[0])
    }

    companion object {
        /** Полиномиальная phi^B(t) = (1, t, t^2)^T. */
        val B = GeneratingSystem(
            name = "B",
            rho = { t -> t }, sigma = { t -> t * t },
            rhoD = { 1.0 }, sigmaD = { t -> 2.0 * t },
            rhoDD = { 0.0 }, sigmaDD = { 2.0 },
        )

        /** Гиперболическая phi^H(t) = (1, sinh t, cosh t)^T. */
        val H = GeneratingSystem(
            name = "H",
            rho = { t -> Math.sinh(t) }, sigma = { t -> Math.cosh(t) },
            rhoD = { t -> Math.cosh(t) }, sigmaD = { t -> Math.sinh(t) },
            rhoDD = { t -> Math.sinh(t) }, sigmaDD = { t -> Math.cosh(t) },
        )

        /** Тригонометрическая phi^T(t) = (1, sin t, cos t)^T. */
        val T = GeneratingSystem(
            name = "T",
            rho = { t -> Math.sin(t) }, sigma = { t -> Math.cos(t) },
            rhoD = { t -> Math.cos(t) }, sigmaD = { t -> -Math.sin(t) },
            rhoDD = { t -> -Math.sin(t) }, sigmaDD = { t -> -Math.cos(t) },
        )
    }
}
