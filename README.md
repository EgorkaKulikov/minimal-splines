# minimal-splines

A library of quadratic minimal splines and quasi-interpolation for the JVM platform. Splines are built from a generating
system φ = (1, ρ, σ): polynomial (1, t, t²), hyperbolic (1, sinh t, cosh t), trigonometric (1, sin t, cos t)
or user-defined. The basis ω_j is determined by the approximation relation Σ_j a_j ω_j(t) = φ(t); every function
ω_j is supported on three adjacent grid intervals, the basis forms a partition of unity and reproduces span φ exactly.

An approximation of a function g is built as the quasi-interpolant P_χ g = Σ_j χ_j(g) ω_j with local approximation
functionals χ_j. Five families are implemented: projection functionals θ, de Boor–Fix functionals ξ^⟨r⟩, their discrete
counterparts ξ̃, averaging functionals μ and three-point functionals λ. Every coefficient χ_j(g) depends on g only on the
support [x_j, x_{j+3}]; the families θ and ξ are biorthogonal to the basis and therefore define projectors onto the spline
space.

## Rationale

Widely used spline libraries implement the polynomial case only. SciPy (`interpolate.BSpline`, `make_lsq_spline`),
Hipparchus (`PolynomialSplineFunction`, `SplineInterpolator`) and Apache Commons Math build piecewise polynomial functions;
generating systems with other ρ and σ, reproducing hyperbolic or trigonometric functions exactly, are not available
there. For the system (1, t, t²) minimal splines coincide with quadratic B-splines, and the library
reproduces them to rounding accuracy.

The interpolating splines of those libraries require solving a global system of linear algebraic equations over
all grid nodes. Quasi-interpolation by local functionals avoids it: the coefficient of ω_j is computed from the
values of g (and, for the family ξ, of its derivative) at a few points of the support, so a change of g on a single interval
affects only the neighbouring coefficients, while the convergence order O(h³) on smooth functions is preserved.

The methods for constructing the basis and the functionals are taken from the works of Yu. K. Demyanovich on minimal splines
and from the works of A. A. Makarov and E. K. Kulikov; the list of sources linked to the code elements is given in `docs/REFERENCES.md`.

## Usage

JDK 21 or newer is required. The dependency `io.github.egorkakulikov:minimal-splines:1.1.0` is published to GitHub Packages;
linear algebra is performed by the `io.github.egorkakulikov:numerical-core` library, which is resolved transitively from
its own repository. Reading both repositories requires a token with the `read:packages` scope (`gpr.user` and `gpr.token`
in `~/.gradle/gradle.properties`, or the corresponding environment variables):

```kotlin
maven {
    url = uri("https://maven.pkg.github.com/EgorkaKulikov/minimal-splines")
    credentials {
        username = providers.gradleProperty("gpr.user").orNull ?: System.getenv("GITHUB_ACTOR")
        password = providers.gradleProperty("gpr.token").orNull ?: System.getenv("GITHUB_TOKEN")
    }
}
maven {
    url = uri("https://maven.pkg.github.com/EgorkaKulikov/numerical-core")
    credentials {
        username = providers.gradleProperty("gpr.user").orNull ?: System.getenv("GITHUB_ACTOR")
        password = providers.gradleProperty("gpr.token").orNull ?: System.getenv("GITHUB_TOKEN")
    }
}
```

## Documentation

The problem statement, evaluation in local interval coordinates, the measured convergence orders and the limits of
applicability are described in `docs/ACCURACY.md`; the sources of the methods are listed in `docs/REFERENCES.md`; the API
documentation is built with `./gradlew dokkaHtml`.

## License

The library is distributed under the Apache License 2.0; the copyright holder is Egor Konstantinovich Kulikov.
The text of the license is in the `LICENSE` file, third-party component notices are in the `NOTICE` file.
