# Базовая линия производительности minimal-splines 0.1.0

Снято до переработки алгоритмов (ветка `rework/1.0`, исходная точка 72f3419). Измерено: построение `MinimalSplineBasis(sys, grid)` (basis),
`projectorCoeffs(f, fD, fDD)` для `ProjFunctionals` (theta), `DeBoorFixFunctionals(basis, 1)` (xi1), `AveragingFunctionals` (mu)
и 10⁴ вызовов `evalSpline` на `java.util.Random(42)` (eval×10⁴); `f(t) = exp(sin 3t)`, сетка `Grid.uniform(n, 0, 1)`; медиана из 5 после 3 прогревов, мс.
Машина: Apple M1 Pro (8 процессоров), Corretto 21.0.10; BLAS/LAPACK — netlib JNILAPACK (системная нативная реализация, Accelerate).
Повторить: `./gradlew benchmark --offline -Pbench.args="100 1000 3000"`.

Отклонение от плана: n = 10⁴ и n = 5000 в 0.1.0 не строятся — `MinimalSplineBasis` падает в `invert3`
(`matrix is singular or near-singular`, det/scale < 1e-12 при h = 1e-4; при h = 2e-4 порог тоже не пройден); максимальный работающий размер из ряда — 3000.

backend: netlib JNILAPACK (нативная BLAS/LAPACK системы)
processors: 8
jvm: OpenJDK 64-Bit Server VM 21.0.10
date: 2026-09-09T19:57:27.601531+03:00

| система | n | basis | theta | xi1 | mu | eval×10⁴ |
|---|---|---|---|---|---|---|
| B | 100 | 0.17 | 0.55 | 0.12 | 0.25 | 0.78 |
| H | 100 | 0.09 | 0.31 | 0.06 | 0.23 | 1.04 |
| B | 1000 | 0.26 | 2.36 | 0.27 | 1.64 | 0.67 |
| H | 1000 | 0.23 | 2.41 | 0.16 | 1.69 | 0.77 |
| B | 3000 | 0.11 | 6.52 | 0.20 | 4.32 | 0.76 |
| H | 3000 | 0.33 | 7.25 | 0.33 | 4.71 | 0.84 |
