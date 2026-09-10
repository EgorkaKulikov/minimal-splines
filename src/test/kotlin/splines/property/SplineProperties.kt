package splines.property

import net.jqwik.api.Arbitrary
import net.jqwik.api.ForAll
import net.jqwik.api.Property
import net.jqwik.api.Provide
import net.jqwik.api.Tag
import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import splines.functionals.AveragingFunctionals
import splines.functionals.DeBoorFixFunctionals
import splines.functionals.DiscreteDeBoorFixFunctionals
import splines.functionals.FunctionalFamily
import splines.functionals.ProjFunctionals
import splines.functionals.ThreePointFunctionals
import splines.metrics.errorEh
import splines.testFunction
import kotlin.math.abs
import kotlin.math.max
import kotlin.random.Random
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Свойства квадратичных минимальных сплайнов на случайных сетках (jqwik): разбиение единицы,
 * точность квазипроекторов на span φ (θ, ξ, μ, λ), биортогональность θ и ξ, идемпотентность,
 * инвариантность базиса к сдвигу (и масштабу для B), носитель и знак ω_j.
 */
@Tag("fast")
class SplineProperties {

    @Provide fun gridCase(): Arbitrary<GridCase> = Generators.gridCases()
    @Provide fun coeffs(): Arbitrary<DoubleArray> = Generators.coefficients()

    private fun samplePoints(grid: Grid, count: Int, seed: Int): DoubleArray {
        val rnd = Random(seed)
        return DoubleArray(count) { grid.a + (grid.b - grid.a) * rnd.nextDouble() }
    }

    private fun families(basis: MinimalSplineBasis): List<FunctionalFamily> = listOf(
        ProjFunctionals(basis),
        DeBoorFixFunctionals(basis, 1),
        DiscreteDeBoorFixFunctionals(basis, 1),
        AveragingFunctionals(basis),
        ThreePointFunctionals(basis),
    )

    @Property(tries = 50)
    fun partitionOfUnity(@ForAll("gridCase") case: GridCase) {
        val basis = MinimalSplineBasis(case.sys, case.grid)
        for (t in samplePoints(case.grid, 50, 17)) {
            var sum = 0.0
            for (j in -2 until case.grid.n) sum += basis.omega(j, t)
            assertTrue(abs(sum - 1.0) <= 1e-12, "$case: Σω_j($t) − 1 = ${sum - 1.0}")
        }
    }

    /**
     * Точность на span φ проверяется для четырёх семейств, точных по построению: θ, ξ, μ, λ.
     * Семейство ξ̃ заменяет f'(x_{j+1}) разделённой разностью по узлам x_j, x_{j+2} и на span φ
     * не точно: на неравномерной сетке погрешность коэффициента есть w·f''·(x_j + x_{j+2} − 2x_{j+1})/2,
     * то есть O(h²) даже для B; для ξ̃ проверяется порядок сходимости (ConvergenceOrderTest).
     */
    @Property(tries = 50)
    fun exactOnSpanPhi(@ForAll("gridCase") case: GridCase, @ForAll("coeffs") c: DoubleArray) {
        val sys = case.sys
        val grid = case.grid
        val basis = MinimalSplineBasis(sys, grid)
        val f: (Double) -> Double = { t -> c[0] + c[1] * sys.rho(t) + c[2] * sys.sigma(t) }
        val fD: (Double) -> Double = { t -> c[1] * sys.rhoD(t) + c[2] * sys.sigmaD(t) }
        val fDD: (Double) -> Double = { t -> c[1] * sys.rhoDD(t) + c[2] * sys.sigmaDD(t) }
        val scale = max(1.0, samplePoints(grid, 200, 3).maxOf { abs(f(it)) })
        for (family in families(basis).filter { it !is DiscreteDeBoorFixFunctionals }) {
            val gD: (Double) -> Double = if (family.usesDerivative) fD else { _ -> 0.0 }
            val gDD: (Double) -> Double = if (family.usesSecondDerivative) fDD else { _ -> 0.0 }
            val coef = family.projectorCoeffs(f, gD, gDD)
            val err = errorEh(f, { t -> basis.evalSpline(coef, t) }, grid)
            assertTrue(err <= 1e-9 * scale, "$case, ${family.name}: E_h = $err, scale = $scale, c = ${c.toList()}")
        }
    }

    @Property(tries = 50)
    fun biorthogonality(@ForAll("gridCase") case: GridCase) {
        val basis = MinimalSplineBasis(case.sys, case.grid)
        val n = case.grid.n
        for (family in families(basis).filter { it.isProjector }) {
            for (i in -2 until n) for (j in -2 until n) {
                val v = family.chi(i).apply(
                    { t -> basis.omega(j, t) },
                    { t -> basis.omegaDeriv(j, t) },
                    { t -> basis.omegaDeriv2(j, t) },
                )
                val expected = if (i == j) 1.0 else 0.0
                assertTrue(abs(v - expected) <= 1e-10, "$case, ${family.name}: χ_$i(ω_$j) = $v")
            }
        }
    }

    @Property(tries = 50)
    fun idempotency(@ForAll("gridCase") case: GridCase) {
        val basis = MinimalSplineBasis(case.sys, case.grid)
        val tf = testFunction()
        for (family in families(basis).filter { it.isProjector }) {
            val gD: (Double) -> Double = if (family.usesDerivative) tf.fD else { _ -> 0.0 }
            val gDD: (Double) -> Double = if (family.usesSecondDerivative) tf.fDD else { _ -> 0.0 }
            val c1 = family.projectorCoeffs(tf.f, gD, gDD)
            val c2 = family.projectorCoeffs(
                { t -> basis.evalSpline(c1, t) },
                { t -> basis.evalSplineDeriv(c1, t) },
                { t -> basis.evalSplineDeriv2(c1, t) },
            )
            val scale = max(1.0, c1.maxOf { abs(it) })
            val diff = c1.indices.maxOf { abs(c1[it] - c2[it]) }
            assertTrue(diff <= 1e-10 * scale, "$case, ${family.name}: max|c₁ − c₂| = $diff, scale = $scale")
        }
    }

    @Property(tries = 50)
    fun shiftInvariance(
        @ForAll("gridCase") case: GridCase,
        @ForAll("shifts") s: Double,
        @ForAll("scales") lambda: Double,
    ) {
        val grid = case.grid
        val n = grid.n
        val basis = MinimalSplineBasis(case.sys, grid)
        val pts = samplePoints(grid, 50, 29)
        val shifted = MinimalSplineBasis(case.sys, Grid(n, DoubleArray(n + 1) { grid.x(it) + s }))
        // Аргумент t + s и узлы x_k + s округляются независимо: разность u = (t+s) − (x_k+s)
        // отличается от t − x_k на величину ≈ ε·(|s| + |t|), а ω_j' = O(1/h_min). Это ошибка
        // представления входа, а не свойство метода, поэтому допуск учитывает её явно.
        val hMin = (0 until n).minOf { grid.x(it + 1) - grid.x(it) }
        val tolShift = 1e-12 + 16.0 * Math.ulp(1.0) * (abs(s) + abs(grid.a) + abs(grid.b)) / hMin
        for (j in -2 until n) for (t in pts) {
            val v0 = basis.omega(j, t)
            val v1 = shifted.omega(j, t + s)
            assertTrue(abs(v0 - v1) <= tolShift, "$case, сдвиг s=$s: ω_$j($t)=$v0, ω_$j(t+s)=$v1, допуск $tolShift")
        }
        if (case.sys === GeneratingSystem.B) {
            val scaled = MinimalSplineBasis(case.sys, Grid(n, DoubleArray(n + 1) { grid.x(it) * lambda }))
            for (j in -2 until n) for (t in pts) {
                val v0 = basis.omega(j, t)
                val v1 = scaled.omega(j, t * lambda)
                // Округление λ·t и λ·x_k независимо: та же ошибка представления входа, что и при сдвиге.
                val tolScale = 1e-12 + 16.0 * Math.ulp(1.0) * (abs(grid.a) + abs(grid.b)) / hMin
                assertTrue(abs(v0 - v1) <= tolScale, "$case, масштаб λ=$lambda: ω_$j($t)=$v0, ω_$j(λt)=$v1, допуск $tolScale")
            }
        }
    }

    @Provide fun shifts(): Arbitrary<Double> = net.jqwik.api.Arbitraries.doubles().between(-100.0, 100.0)
    @Provide fun scales(): Arbitrary<Double> = net.jqwik.api.Arbitraries.doubles().between(1e-3, 1e3).ofScale(3)

    @Property(tries = 50)
    fun supportAndSign(@ForAll("gridCase") case: GridCase) {
        val grid = case.grid
        val n = grid.n
        val basis = MinimalSplineBasis(case.sys, grid)
        val rnd = Random(41)
        for (j in -2 until n) {
            val lo = grid.x(j)
            val hi = grid.x(j + 3)
            // Вне носителя (внутри [a, b]) — ровно 0.0.
            var outside = 0
            var attempts = 0
            while (outside < 20 && attempts < 2000) {
                attempts++
                val t = grid.a + (grid.b - grid.a) * rnd.nextDouble()
                if (t >= lo && t <= hi) continue
                outside++
                assertEquals(0.0, basis.omega(j, t), "$case: ω_$j($t) ≠ 0 вне [$lo, $hi]")
            }
            // Внутри носителя — неотрицательность; для B — строгая положительность.
            for (q in 1 until 40) {
                val t = lo + (hi - lo) * q / 40.0
                if (t < grid.a || t > grid.b || t <= lo || t >= hi) continue
                val v = basis.omega(j, t)
                assertTrue(v >= -1e-14, "$case: ω_$j($t) = $v < 0 внутри [$lo, $hi]")
                if (case.sys === GeneratingSystem.B) assertTrue(v > 0.0, "$case: ω_$j($t) = $v не > 0")
            }
        }
    }
}
