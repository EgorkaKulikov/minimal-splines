package splines.golden

import numerics.backend.Backends
import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import splines.functionals.AveragingFunctionals
import splines.functionals.DeBoorFixFunctionals
import splines.functionals.DiscreteDeBoorFixFunctionals
import splines.functionals.FunctionalFamily
import splines.functionals.ProjFunctionals
import splines.functionals.ThreePointFunctionals
import splines.golden.GoldenIo.dbl
import splines.golden.GoldenIo.mat
import splines.golden.GoldenIo.vec
import splines.metrics.errorEh
import splines.nonDegenerate

/**
 * Computation of the reference quantities. The same code produces the golden reference (GoldenGenerate)
 * and repeats the computation during verification (Golden*Test): trees of values are compared.
 * Numbers are encoded as hexadecimal raw bits (see GoldenIo).
 */
object GoldenCompute {
    fun meta(): Map<String, Any?> = linkedMapOf(
        "version" to (System.getProperty("golden.version") ?: "0.1.0"),
        "backend" to Backends.default().name,
        "date" to java.time.Instant.now().toString(),
    )

    // ---- grid.json -------------------------------------------------------------

    fun gridCase(g: Grid): Map<String, Any?> = linkedMapOf(
        "n" to g.n,
        "x" to vec(DoubleArray(g.n + 5) { g.x(it - 2) }),
        "h" to dbl(g.h),
        "breakpoints" to vec(g.breakpoints),
        "isCoincident" to (-2..g.n + 1).map { g.isCoincident(it) },
        "breakpointInclusionEps" to dbl(g.breakpointInclusionEps),
    )

    // ---- basis.json ------------------------------------------------------------

    /** Derivatives of the basis functions are written only for n <= 8 (to save space). */
    fun basisCase(sys: GeneratingSystem, g: Grid): Map<String, Any?> {
        val basis = MinimalSplineBasis(sys, g)
        val t = GoldenInputs.controlPoints(g)
        val c = GoldenInputs.coeffs(g.n)
        val js = -2..g.n - 1
        fun rows(f: (Int, Double) -> Double): List<List<String>> =
            mat(js.map { j -> DoubleArray(t.size) { f(j, t[it]) } }.toTypedArray())
        val m = linkedMapOf<String, Any?>(
            "n" to g.n,
            "t" to vec(t),
            "interval" to t.map { basis.interval(it) },
            "omega" to rows { j, x -> basis.omega(j, x) },
            "evalSpline" to vec(DoubleArray(t.size) { basis.evalSpline(c, t[it]) }),
            "evalSplineDeriv" to vec(DoubleArray(t.size) { basis.evalSplineDeriv(c, t[it]) }),
            "evalSplineDeriv2" to vec(DoubleArray(t.size) { basis.evalSplineDeriv2(c, t[it]) }),
        )
        if (g.n <= 8) {
            m["omegaDeriv"] = rows { j, x -> basis.omegaDeriv(j, x) }
            m["omegaDeriv2"] = rows { j, x -> basis.omegaDeriv2(j, x) }
        }
        return m
    }

    // ---- functionals.json ------------------------------------------------------

    val familyKeys: List<String> = listOf(
        "proj", "deBoor-0", "deBoor-1", "deBoor-2",
        "discreteDeBoor-1", "discreteDeBoor-2", "averaging", "threePoint",
    )

    fun family(key: String, basis: MinimalSplineBasis): FunctionalFamily = when (key) {
        "proj" -> ProjFunctionals(basis)
        "deBoor-0" -> DeBoorFixFunctionals(basis, 0)
        "deBoor-1" -> DeBoorFixFunctionals(basis, 1)
        "deBoor-2" -> DeBoorFixFunctionals(basis, 2)
        "discreteDeBoor-1" -> DiscreteDeBoorFixFunctionals(basis, 1)
        "discreteDeBoor-2" -> DiscreteDeBoorFixFunctionals(basis, 2)
        "averaging" -> AveragingFunctionals(basis)
        "threePoint" -> ThreePointFunctionals(basis)
        else -> error("unknown family: $key")
    }

    /**
     * If building the family or the computation throws, its exception class is recorded:
     * that is part of the 0.1.0 behaviour as well (for example, families on too short grids).
     */
    fun familyCase(key: String, basis: MinimalSplineBasis): Map<String, Any?> {
        val m = linkedMapOf<String, Any?>()
        try {
            val fam = family(key, basis)
            m["isProjector"] = fam.isProjector
            m["usesDerivative"] = fam.usesDerivative
            m["usesSecondDerivative"] = fam.usesSecondDerivative
            m["cChi"] = dbl(fam.cChi())
            m["coeffs"] = GoldenInputs.functions.associate { it.name to vec(fam.projectorCoeffs(it.f, it.fD, it.fDD)) }
            if (fam is ProjFunctionals && basis.n >= 3) {
                try {
                    val f1 = GoldenInputs.functions[0]
                    val js = 0..basis.n - 3
                    m["closedFormAbsSum"] = vec(js.map { fam.closedFormInternal(it).absSum() }.toDoubleArray())
                    m["closedFormApplyF1"] = vec(js.map { fam.closedFormInternal(it).apply(f1.f, f1.fD, f1.fDD) }.toDoubleArray())
                } catch (e: Exception) {
                    m["closedFormError"] = e::class.simpleName
                }
            }
        } catch (e: Exception) {
            return linkedMapOf("error" to e::class.simpleName)
        }
        return m
    }

    fun functionalsCase(sys: GeneratingSystem, g: Grid): Map<String, Any?> {
        val basis = MinimalSplineBasis(sys, g)
        return familyKeys.associateWith { familyCase(it, basis) }
    }

    // ---- metrics.json ----------------------------------------------------------

    val metricsGridKeys: List<String> = listOf("uniform", "graded").flatMap { k -> listOf(4, 8, 16, 32).map { "$k-$it" } }

    fun errorEhCases(): Map<String, Any?> = linkedMapOf<String, Any?>().apply {
        for ((sName, sys) in GoldenInputs.systems) for (gKey in metricsGridKeys) {
            val g = GoldenInputs.grids.getValue(gKey)
            val basis = MinimalSplineBasis(sys, g)
            val theta = ProjFunctionals(basis)
            for (fn in GoldenInputs.functions) {
                val c = theta.projectorCoeffs(fn.f, fn.fD, fn.fDD)
                put("$sName-$gKey-${fn.name}", dbl(errorEh(fn.f, { t -> basis.evalSpline(c, t) }, g)))
            }
        }
    }

    // ---- algebra.json ----------------------------------------------------------

    fun nonDegenerateCases(): Map<String, Any?> =
        GoldenInputs.grids.mapValues { (_, g) -> (-2..g.n - 1).map { nonDegenerate(g, it) } }

    fun phiCases(): Map<String, Any?> = GoldenInputs.systems.mapValues { (_, s) ->
        val p = GoldenInputs.phiPoints
        linkedMapOf(
            "phi" to mat(p.map { s.phi(it) }.toTypedArray()),
            "phiD" to mat(p.map { s.phiD(it) }.toTypedArray()),
            "phiDD" to mat(p.map { s.phiDD(it) }.toTypedArray()),
            "wronskian" to vec(DoubleArray(p.size) { s.wronskian(p[it]) }),
        )
    }
}
