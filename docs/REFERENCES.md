# Provenance of the methods

This document records the provenance of every method of the library: of the mathematical construction,
of the formula or of the engineering decision. Each table lists the element, the symbol that implements it
and the status of the construction.

## Status conventions

- **Classical** — the construction is known from the literature; the implementation is checked against
  the published formulas and against the tests.
- **Original work** — an engineering decision of the library with no external source; the
  justification is given in the KDoc of the corresponding element and in `docs/ACCURACY.md`.

## 1. Generating systems and basis

| Element | Implementation | Status |
|---|---|---|
| Generating system φ = (1, ρ, σ) with a non-degenerate Wronskian; polynomial system (1, t, t²) | `GeneratingSystem`, `GeneratingSystem.B` | Classical |
| Hyperbolic system (1, sinh t, cosh t) | `GeneratingSystem.H` | Classical |
| Trigonometric system (1, sin t, cos t) | `GeneratingSystem.T` | Classical |
| Grid with triple nodes at the ends of the interval | `Grid` | Classical |
| Approximation relation Σ_j a_j ω_j = φ, vectors a_j from the values of φ, φ' at the nodes, basis ω(t) = M_k⁻¹ φ(t) | `MinimalSplineBasis`, `MinimalSplineBasis.computeA` | Classical |
| Explicit formulas for the basis functions of the systems B and H (oracle in the tests) | `src/test/kotlin/splines/ReferenceSplines.kt` | Classical |
| Coincidence with quadratic B-splines for the system B; Cox–de Boor recursion in the tests | `src/test/kotlin/splines/DeBoorOracleTest.kt` | Classical |
| Local interval coordinates ψ_k = T_k φ; covariance a_j^{Tφ} = T a_j^{φ} | `LocalFrame`, `GeneratingSystem.localFrame`, `MinimalSplineBasis.localApproximationMatrix` | Original work |
| Degeneracy criterion based on the condition number of T_k M_k, and the relative significance threshold for the denominators | `MinimalSplineBasis.MAX_CONDITION`, `DEGENERACY_RELATIVE_EPS`, `nonDegenerate` | Original work |
| Tangent vector of the local frame: the vectors a_j are built from it and do not depend on a scalar factor of φ' | `LocalFrame.tangent`, `MinimalSplineBasis.computeA` | Original work |
| Reparametrized system φ = (1, g, g²): ω_j = B_j ∘ g, where B_j are the quadratic B-splines on the knots g(x_j); established symbolically and verified by `ReparametrizedSystemTest` | `Reparametrization`, `Reparametrization.power`, `GeneratingSystem.reparametrized` | Original work |
| Power-graded grid x_j = a + (b − a)(j/n)^r (H. Brunner, The numerical solution of weakly singular Volterra integral equations by collocation on graded meshes, Math. Comp. 45 (1985), no. 172, 417–437, doi:10.1090/S0025-5718-1985-0804933-3); the symmetric variant applies the same law on each half of the interval | `Grid.power`, `Grid.symmetricPower` | Classical |

## 2. Approximation functionals

| Element | Implementation | Status |
|---|---|---|
| Projection functionals θ_j: five-point closed form on the nodes x_j, x_{j+3} and the midpoints of the support intervals; biorthogonality θ_i(ω_j) = δ_ij | `ProjFunctionals` | Classical |
| De Boor–Fix functionals ξ^⟨r⟩, r ∈ {0, 1, 2}, for polynomial splines | `DeBoorFixFunctionals` | Classical |
| Functionals ξ^⟨r⟩ for minimal splines built from an arbitrary generating system; biorthogonality ξ_i(ω_j) = δ_ij | `DeBoorFixFunctionals` | Classical |
| Boundary functionals θ, ξ, ξ̃ at j = −2 and j = n − 1 as the values g(x₀), g(x_n) | `ProjFunctionals`, `DeBoorFixFunctionals`, `DiscreteDeBoorFixFunctionals` | Original work |
| Discrete de Boor–Fix functionals ξ̃^⟨r⟩, r ∈ {1, 2}: the derivative is replaced by a divided difference over the neighbouring nodes | `DiscreteDeBoorFixFunctionals` | Original work |
| Averaging functionals μ_j from the values at the points of the auxiliary grid y_j = x_{j+1} + θ (x_{j+2} − x_{j+1}) | `AveragingFunctionals` | Classical |
| Three-point functionals λ_j from the values at x_{j+1}, x_{j+3/2}, x_{j+2} | `ThreePointFunctionals` | Classical |
| Evaluation of the coefficients of the functionals through the local representation ψ_k instead of φ | `FunctionalFamily` and its subclasses | Original work |
| Placement of the midpoints of θ_j by a sampling rule: arithmetic (l + r)/2, or g⁻¹((g(l) + g(r))/2) for a reparametrized system | `ThetaSampling`, `ThetaSampling.ARITHMETIC`, `ThetaSampling.reparametrized`, `ProjFunctionals` | Original work |

## 3. Metrics

| Element | Implementation | Status |
|---|---|---|
| Error E_h as the maximum of \|g − P_χ g\| over a uniform control grid with a given number of points per interval | `errorEh`, `DEFAULT_CONTROL_REFINEMENT` | Classical |
| Observed convergence order log₂(E_h / E_{h/2}) | convergence tests; the functions `orders`, `reliableOrders` of the numerical-core library | Classical |
