# Rules for changing minimal-splines

minimal-splines is a library of quadratic minimal splines and quasi-interpolation written in
Kotlin/JVM. The rules below govern every change to the code and the documentation.

## Language

KDoc, comments, exception messages, commit messages and all documentation are written in English;
the only exception is `docs/ABSTRACT.md`, the software-registration abstract, which is kept in
Russian. Identifiers in the code are English. The style of the text is descriptive: no words in
capitals, no jargon, no references to external projects or to the version history.

## Public API

- `explicitApi()` is enabled: every new element is declared `internal`. An element becomes public
  only once it has KDoc and a test covering its contract.
- numerical-core types (`NumericsContext`, `DenseMatrix`) appear in the signatures of the public API,
  so the dependency is declared in the `api` scope; no other external dependencies are added to the
  main source set.
- The behaviour of a method (checks, exceptions, shape of the result) does not depend on the
  BLAS/LAPACK implementation.

## Numerical methods

- Any quantity of the generating system on a grid interval is computed only through the local
  representation `basis.frame(k).psi`, `psiD`, `psiDD`. Accessing `sys.rho`, `sys.sigma` and their
  derivatives in the global variable t inside the basis and functional code is not allowed: it
  causes a loss of significance that depends on the length and the position of the interval (the
  justification is in `docs/ACCURACY.md`, section "Local interval coordinates").
- Every family of functionals comes with tests of biorthogonality or of exactness on
  `span{1, rho, sigma}`, with a comparison against the closed formula of the source and with an
  entry in `docs/REFERENCES.md`.
- The thresholds `MinimalSplineBasis.MAX_CONDITION` and `DEGENERACY_RELATIVE_EPS` are changed only
  with a justification in `docs/ACCURACY.md`.
- The golden references in `src/test/resources/golden` are regenerated with
  `./gradlew regenerateGolden` only when the behaviour is changed intentionally, with the reason
  recorded in `src/test/resources/golden/README.md`.
- A counterexample found by the jqwik properties is eliminated by fixing the code; relaxing a
  tolerance or narrowing a generator is not allowed.

## Commands

    ./gradlew build                                  # tests, Kover coverage threshold
    ./gradlew test -Dnumerics.backend=java           # portable Java implementation
    ./gradlew test -Dnumerics.backend=native         # system BLAS/LAPACK implementation
    ./gradlew benchmark -Pbench.args="100 1000"      # performance measurements
    ./gradlew dokkaHtml                              # API documentation

The `check` task includes the Kover coverage threshold check; a change that drops the coverage below
the threshold fails the build. All commands are run with the `--offline` flag.

Releasing a version: the tag `vX.Y.Z` triggers publication to GitHub Packages from CI.

## Documentation

After an API change, `README.md` and the corresponding document in `docs/` are updated. The numbers in
`docs/ACCURACY.md` are taken from the reports `build/reports/*.tsv` produced by the `test` task.

## Commits

Commit messages are written in English with an area prefix: `api:`, `algo:`, `test:`, `docs:`,
`build:`. One change corresponds to one commit; the message describes what changed and why.
