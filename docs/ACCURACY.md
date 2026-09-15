# Accuracy and limits of applicability

This document describes the mathematical problem statement, the evaluation scheme that makes the
conditioning independent of the grid and of the interval, the measured convergence orders of the
quasi-interpolants, and the limits beyond which the construction terminates with an exception. All the
numbers come from the tests of the library; the corresponding tables are produced by `./gradlew test`
in the `build/reports` directory (the `*.tsv` files named below in each section).

## Problem statement

A generating system is a vector function φ(t) = (1, ρ(t), σ(t))ᵀ with a non-degenerate Wronskian
det(φ, φ', φ''). The built-in systems are `GeneratingSystem.B` (1, t, t²), `H` (1, sinh t, cosh t) and
`T` (1, sin t, cos t); a user-defined system is given by the functions ρ, σ and their derivatives up to
the second order. A `Grid` on the interval [a, b] consists of n intervals with triple nodes at the ends:
x₋₂ = x₋₁ = x₀ = a, x_n = x_{n+1} = x_{n+2} = b.

The basis of quadratic minimal splines {ω_j}, j = −2, …, n − 1, is determined by the approximation
relation

    Σ_j a_j ω_j(t) = φ(t),   t ∈ [a, b],

in which the vectors a_j ∈ ℝ³ are built from the values of φ and φ' at the nodes x_{j+1}, x_{j+2}: a_j is
the direction vector of the intersection of the planes span{φ(x_{j+1}), φ'(x_{j+1})} and span{φ(x_{j+2}), φ'(x_{j+2})},
normalized by the condition a_j = φ(x_{j+1}) − c φ'(x_{j+1}). On the interval (x_k, x_{k+1}) exactly three
functions ω_{k−2}, ω_{k−1}, ω_k are non-zero, and the relation takes the form M_k ω(t) = φ(t) with the matrix
M_k = (a_{k−2} | a_{k−1} | a_k), whence ω(t) = M_k⁻¹ φ(t). The support of ω_j is [x_j, x_{j+3}]; the first
component φ₀ ≡ 1 yields the partition of unity Σ_j ω_j ≡ 1; the splines belong to C¹[a, b] and reproduce
span{1, ρ, σ} exactly.

The quasi-interpolant of a function g has the form P_χ g = Σ_j χ_j(g) ω_j, where the χ_j are linear functionals
that depend on g only on the support [x_j, x_{j+3}]. The families θ (`ProjFunctionals`) and ξ^⟨r⟩
(`DeBoorFixFunctionals`, r ∈ {0, 1, 2}) are biorthogonal to the basis, χ_i(ω_j) = δ_ij, and therefore P_χ is
a projector onto the spline space. The families ξ̃ (`DiscreteDeBoorFixFunctionals`), μ
(`AveragingFunctionals`) and λ (`ThreePointFunctionals`) are not projectors, but μ and λ are exact on
span{1, ρ, σ}. The error is measured by the function `errorEh`: the maximum of |g − P_χ g| over a uniform
control grid with `DEFAULT_CONTROL_REFINEMENT = 100` points per interval.

## Local interval coordinates

Inverting M_k directly in the global variable t is unusable for large n and for intervals that are far from
zero or whose length differs from one. The columns of M_k contain the values of φ at neighbouring nodes and, for a small
step h, are nearly collinear: for the system B on [0, 1] the condition number grows as n² — from 1.7·10³
at n = 10 to 1.8·10⁹ at n = 10⁴, and on the intervals [0, 10⁶] and [0, 10⁻⁶] it reaches 3·10¹⁶, that is,
the level at which the inverse matrix contains no correct digits. The functionals ξ, computed from differences
of the values of ρ, σ, lost biorthogonality as 1/L² with the interval length L: for the system H on [0, 10⁻³] the
biorthogonality residual was 3·10⁻⁸, and on [100, 101] the construction terminated with an exception
because of the mutual cancellation of terms of order e¹⁰⁰ (files `conditioning-diagnostic.tsv`,
`conditioning-diagnostic-segments.tsv`, `xi-diagnostic.tsv`).

The cause lies not in the problem itself but in the choice of coordinates. The approximation relation is invariant
under the substitution φ → Tφ with a non-singular matrix T: multiplying from the left gives (T M_k) ω(t) = Tφ(t) with
the same solution ω(t), and the vector a_j is covariant, a_j^{Tφ} = T a_j^{φ}, because the coefficient c is
the ratio of two scalar products with one and the same normal φ(x_{j+2}) × φ'(x_{j+2}). Therefore, on
every interval (x_k, x_{k+1}) the local representation ψ_k = T_k φ is used
(`GeneratingSystem.localFrame`, class `LocalFrame`), whose components are of order one:

| System | ψ_k(t), u = t − x_k, h = x_{k+1} − x_k, l = min(h, 1) |
|---|---|
| B | (1, s, s²), s = u/h |
| H | (1, sinh(u)/l, 2 sinh²(u/2)/l²) |
| T | (1, sin(u)/l, 2 sin²(u/2)/l²) |

The third components for H and T are written through the square of the half argument rather than as (cosh u − 1)/l²
and (1 − cos u)/l², in order to avoid the subtraction of close numbers. The matrix T_k M_k = (T_k a_{k−2} | T_k a_{k−1} | T_k a_k)
is assembled by the same formula for a_j, written for ψ_k instead of φ, and is inverted by numerical-core.
Its condition number equals 13 for the system B and 21 for the systems H and T for any number of intervals and
any interval on a uniform grid (`conditioning-diagnostic.tsv`). The biorthogonality residual of the functionals ξ in local
coordinates does not exceed 5·10⁻¹⁶ on intervals of length from 10⁻⁴ to 1 (`xi-diagnostic.tsv`), and the relative
error of reproducing a function from span φ by the families θ, ξ, μ and λ at n = 10⁴ on [0, 1] and on the
interval [100, 101] does not exceed 10⁻¹⁵ (`functionals-diagnostic.tsv`).

The degeneracy criterion is an estimate of the condition number of T_k M_k: if it exceeds
`MinimalSplineBasis.MAX_CONDITION = 10⁸`, an exception is raised indicating the interval. The threshold means
that the inverse matrix retains at least eight significant digits out of sixteen; for the built-in
systems on a non-degenerate grid the value 13 or 21 is far from it, and hitting the threshold points to
a degeneracy of the grid itself or of the user-defined system. The representability of ψ_k in
double precision is checked separately: for the system H the quantity sinh(u/l) overflows when u/l > 710, that is, for a step
h > 355, and the exception reports the interval and the node at which the value is non-finite.

A user-defined `GeneratingSystem` for which no local representation is given is evaluated in
global coordinates (T_k = I). For the systems (1, eᵗ, e²ᵗ) on [0, 1] and (1, t, t³) on [1, 2] at n = 50 the
partition of unity holds to an accuracy of 2·10⁻¹² (`boundary-cases.tsv`); the behaviour for large n and on
distant intervals for such systems is determined by the conditioning of M_k in global coordinates.

## Convergence orders

The observed order p = log₂(E_h / E_{h/2}) was measured on the function f(t) = exp(sin 3t), t ∈ [0, 1], at
n = 8, 16, …, 128 on the uniform grid `Grid.uniform` and on the grid `Grid.quasiUniform` with perturbed
nodes; the table gives the values for the pair n = 64 → 128 (`convergence-orders.tsv`). The results
for the systems B, H and T agree to the third digit: on an interval of length h a smooth function
is approximated by an element of each of the spaces span φ with error O(h³), and the difference between the
systems shows up only in the constant.

| Family | Order for the values | Order for the derivative | E_h at n = 128 |
|---|---|---|---|
| θ | 3.003–3.007 | 1.98–2.00 | 4.2·10⁻⁷ |
| ξ = ξ^⟨1⟩ | 2.999–3.004 | 1.98–2.00 | 2.6·10⁻⁶ |
| μ | 3.045–3.081 | — | — |
| λ | 3.018–3.030 | — | — |
| ξ̃, interior intervals | 2.996–2.998 | — | 7·10⁻⁶ |
| ξ̃, whole interval | 1.993–1.999 | ≈ 1.0 | 4.6·10⁻⁵ |

The families θ, ξ, μ and λ have the third order for the values and the second for the first derivative, as follows
from the exactness on span φ for quadratic splines of class C¹. At equal order, the projectors θ and ξ
differ in the constant: the measured error of ξ is six times the error of θ.

## The functionals ξ̃

The family ξ̃^⟨r⟩, r ∈ {1, 2}, is obtained from ξ^⟨r⟩ by replacing the derivative g' with the central divided
difference over the neighbouring grid nodes. The functional uses only the values of g, but is no longer exact
on span φ: on a non-uniform grid the error of the coefficient is of order O(h² |g''|), and on the
interior intervals the observed order remains third.

At the ends of the interval the nodes are multiple, x₋₂ = x₋₁ = x₀, and the central difference degenerates into a one-sided one.
The error of the difference approximation of the derivative becomes a first-order quantity in h, which for
functions from span φ gives a spline deviation of ≈ 0.9975·h² in a boundary layer of width about 3h; over the whole
interval the order drops to the second for the values and to the first for the derivative. This is a property of the
construction, not of the implementation: the measured defect matches the defect predicted by the Taylor expansion of the
one-sided difference. If the third order is required over the whole interval, the families θ, μ or λ, which do not
require the derivative, should be used.

## Limits of applicability

The trigonometric system is applicable on an interval of length less than π. For a greater length, among the points
involved in the construction there are pairs at distance π, for which the determinant built from the derivatives
of ρ and σ, equal to the sine of the difference of the arguments, vanishes; the construction of the basis or of the functionals
terminates with a degeneracy exception.

The hyperbolic system is limited by the double range: for a grid step h > 355 the values of sinh(u/l) are
non-finite, and the construction terminates with an overflow exception indicating the interval. The interval
[0, 10⁶] at n = 100 (h = 10⁴) is not representable for the system H; the same interval for the system B is built
with condition number 13 (`conditioning-diagnostic-segments.tsv`).

The degeneracy of the approximation relation and of the denominators of the functionals is checked by a relative
criterion: a quantity is considered insignificant if it is less than `DEGENERACY_RELATIVE_EPS = 10⁻¹²` of the sum of the
absolute values of the terms it was obtained from. Unlike an absolute threshold, the criterion does not depend on the scale of the
interval, and it fires when nodes coincide outside the ends of the grid and for a degenerate
user-defined system.

## Agreement of the BLAS/LAPACK implementations

The inversion of the matrices T_k M_k and the solution of the small systems of linear algebraic equations in the families
θ, μ and λ are performed by the BLAS/LAPACK implementation selected in `NumericsContext` (the system one or the
portable Java one). The largest relative discrepancy between the implementations is 5.4·10⁻¹⁶ for the
basis and 8.2·10⁻¹⁶ for the functionals (`backend-agreement.tsv`, the systems B, H, T on the uniform and
quasi-uniform grids). The families ξ and ξ̃ are computed by closed formulas without resorting to
linear algebra, and their result does not depend on the implementation.

## Verification of the polynomial case

For the system B, minimal splines coincide with the quadratic B-splines N_{j,2} on the same grid with
triple boundary nodes. The B-splines are computed independently in the tests, by the Cox–de Boor recursion, and are
compared with ω_j on 12 grids (uniform, quasi-uniform, geometric and graded):
the largest discrepancy of the values is 5.6·10⁻¹⁶, of the derivatives 2.8·10⁻¹⁴ (`deboor-oracle.tsv`).
The agreement confirms that the approximation relation, the local coordinates and the matrix inversion
reproduce the classical basis to rounding accuracy.
