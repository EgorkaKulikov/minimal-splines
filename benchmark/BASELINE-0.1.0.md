# Performance baseline of minimal-splines 0.1.0

Taken before the algorithms were reworked (branch `rework/1.0`, starting point 72f3419). Measured: construction of `MinimalSplineBasis(sys, grid)` (basis),
`projectorCoeffs(f, fD, fDD)` for `ProjFunctionals` (theta), `DeBoorFixFunctionals(basis, 1)` (xi1), `AveragingFunctionals` (mu)
and 10⁴ calls of `evalSpline` on `java.util.Random(42)` (eval×10⁴); `f(t) = exp(sin 3t)`, grid `Grid.uniform(n, 0, 1)`; median of 5 after 3 warm-ups, ms.
Machine: Apple M1 Pro (8 processors), Corretto 21.0.10; BLAS/LAPACK — netlib JNILAPACK (the system native implementation, Accelerate).
To reproduce: `./gradlew benchmark --offline -Pbench.args="100 1000 3000"`.

Deviation from the plan: n = 10⁴ and n = 5000 are not built in 0.1.0 — `MinimalSplineBasis` fails in `invert3`
(`matrix is singular or near-singular`, det/scale < 1e-12 at h = 1e-4; at h = 2e-4 the threshold is not met either); the largest working size of the series is 3000.

backend: netlib JNILAPACK (the system native BLAS/LAPACK)
processors: 8
jvm: OpenJDK 64-Bit Server VM 21.0.10
date: 2026-09-09T19:57:27.601531+03:00

| system | n | basis | theta | xi1 | mu | eval×10⁴ |
|---|---|---|---|---|---|---|
| B | 100 | 0.17 | 0.55 | 0.12 | 0.25 | 0.78 |
| H | 100 | 0.09 | 0.31 | 0.06 | 0.23 | 1.04 |
| B | 1000 | 0.26 | 2.36 | 0.27 | 1.64 | 0.67 |
| H | 1000 | 0.23 | 2.41 | 0.16 | 1.69 | 0.77 |
| B | 3000 | 0.11 | 6.52 | 0.20 | 4.32 | 0.76 |
| H | 3000 | 0.33 | 7.25 | 0.33 | 4.71 | 0.84 |
