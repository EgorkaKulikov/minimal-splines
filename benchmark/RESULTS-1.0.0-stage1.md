# Performance of minimal-splines after stage 1 (1c)

Taken on the branch `rework/1.0` after the matrices of the approximation relation and the small linear systems were moved to
`DenseMatrix`/LAPACK of numerical-core 1.0.0 and the degeneracy criterion was changed to an estimate of the condition
number (`MinimalSplineBasis.MAX_CONDITION = 1e8`). The methodology and the machine are as in
`BASELINE-0.1.0.md`: Apple M1 Pro (8 processors), Corretto 21.0.10, netlib JNILAPACK (Accelerate);
grid `Grid.uniform(n, 0, 1)`, median of 5 after 3 warm-ups, ms.
To reproduce: `./gradlew benchmark --offline -Pbench.args="100 1000 2000 3000 10000"`.

## Comparison with 0.1.0

| operation | system | n | 0.1.0 | after 1c | after 1d |
|---|---|---|---|---|---|
| basis | B | 100 | 0.17 | 0.58 | 0.67 |
| basis | H | 100 | 0.09 | 0.42 | 0.48 |
| basis | B | 1000 | 0.26 | 3.48 | 3.69 |
| basis | H | 1000 | 0.23 | 3.06 | 3.44 |
| basis | B | 2000 | — | 5.81 | 5.84 |
| basis | H | 2000 | 0.33 (n = 3000) | not built | 5.99 |
| basis | B, H | 3000 | 0.11 / 0.33 | not built | 8.50 / 8.95 |
| basis | B, H | 10000 | not built (det/scale) | not built (cond > 10⁸) | 28.89 / 29.99 |
| theta | B | 100 | 0.55 | 0.59 | 0.56 |
| theta | H | 100 | 0.31 | 0.50 | 0.66 |
| theta | B | 1000 | 2.36 | 2.59 | 2.52 |
| theta | H | 1000 | 2.41 | 2.37 | 2.38 |
| xi1 | B | 100 | 0.12 | 0.13 | 0.13 |
| xi1 | H | 100 | 0.06 | 0.06 | 0.09 |
| xi1 | B | 1000 | 0.27 | 0.28 | 0.32 |
| xi1 | H | 1000 | 0.16 | 0.38 | 0.33 |
| mu | B | 100 | 0.25 | 0.30 | 0.35 |
| mu | H | 100 | 0.23 | 0.27 | 0.25 |
| mu | B | 1000 | 1.64 | 1.77 | 1.64 |
| mu | H | 1000 | 1.69 | 1.60 | 1.49 |
| eval×10⁴ | B | 100 | 0.78 | 1.06 | 0.76 |
| eval×10⁴ | H | 100 | 1.04 | 0.78 | 0.57 |
| eval×10⁴ | B | 1000 | 0.67 | 0.67 | 0.68 |
| eval×10⁴ | H | 1000 | 0.77 | 0.75 | 0.77 |

Full output of the `benchmark` task:

| system | n | basis | theta | xi1 | mu | eval×10⁴ |
|---|---|---|---|---|---|---|
| B | 100 | 0.58 | 0.59 | 0.13 | 0.30 | 1.06 |
| H | 100 | 0.42 | 0.50 | 0.06 | 0.27 | 0.78 |
| B | 1000 | 3.48 | 2.59 | 0.28 | 1.77 | 0.67 |
| H | 1000 | 3.06 | 2.37 | 0.38 | 1.60 | 0.75 |
| B | 2000 | 5.81 | 4.56 | 0.13 | 2.85 | 0.76 |
| H | 2000 | not built: cond(M_1999) = 1.38·10⁸ | | | | |
| B | 3000 | not built: cond(M_2886) = 1.0007·10⁸ | | | | |
| H | 3000 | not built: cond(M_0) = 1.08·10⁸ | | | | |
| B | 10000 | not built: cond(M_0) = 3.0·10⁸ | | | | |
| H | 10000 | not built: cond(M_0) = 1.2·10⁹ | | | | |

## After stage 1d (local interval coordinates)

Taken after the move to the local representation of the generating system psi_k = T_k phi on every interval
(`GeneratingSystem.localFrame`, `MinimalSplineBasis.localApproximationMatrix`). The same machine, JVM, BLAS/LAPACK implementation and methodology;
`./gradlew benchmark --offline -Pbench.args="100 1000 2000 3000 10000"`.

| system | n | basis | theta | xi1 | mu | eval×10⁴ |
|---|---|---|---|---|---|---|
| B | 100 | 0.67 | 0.56 | 0.13 | 0.35 | 0.76 |
| H | 100 | 0.48 | 0.66 | 0.09 | 0.25 | 0.57 |
| B | 1000 | 3.69 | 2.52 | 0.32 | 1.64 | 0.68 |
| H | 1000 | 3.44 | 2.38 | 0.33 | 1.49 | 0.77 |
| B | 2000 | 5.84 | 4.45 | 0.14 | 2.79 | 0.73 |
| H | 2000 | 5.99 | 4.62 | 0.26 | 2.97 | 0.81 |
| B | 3000 | 8.50 | 6.32 | 0.15 | 3.79 | 0.77 |
| H | 3000 | 8.95 | 6.92 | 0.31 | 4.30 | 0.84 |
| B | 10000 | 28.89 | 21.85 | 0.42 | 12.78 | 0.89 |
| H | 10000 | 29.99 | 24.41 | 0.97 | 14.35 | 0.99 |

Condition number of the matrix being inverted (the largest over k), on the uniform grid [0,1]:

| n | B: cond(M_k) global | B: cond(T_k M_k) | H: cond(M_k) global | H: cond(T_k M_k) | T: cond(T_k M_k) |
|---|---|---|---|---|---|
| 10 | 1.7·10³ | 13.0 | 3.3·10³ | 21.1 | 20.9 |
| 100 | 1.8·10⁵ | 13.0 | 3.4·10⁵ | 21.0 | 21.0 |
| 1000 | 1.8·10⁷ | 13.0 | 3.4·10⁷ | 21.0 | 21.0 |
| 10000 | 1.8·10⁹ | 13.0 | 3.4·10⁹ | 21.0 | 21.0 |

On the intervals [100,101], [0,10⁶], [0,10⁻⁶] at n = 100 the global cond(M_k) for B equals 3·10¹², 3·10¹⁶, 3·10¹⁶,
the local one equals 13.0; for H on [100,101] the local cond equals 21.0 (the global M_k is not computed: the denominator
overflows in global coordinates). The partition of unity in all the listed cases is 2.2·10⁻¹⁶.
The cost of building the basis is linear in n (≈ 3 µs per interval), n = 3·10³ and 10⁴ are built.

## Assessment

1. Building the basis became 3–5 times slower at n = 10² and about 13 times slower at n = 10³: each of the n
   intervals costs two LAPACK calls through JNI (dgecon after the LU factorization, and dgetri),
   so the call overhead dominates over the 3×3 arithmetic. The cost is linear in n —
   about 3 µs per interval (3.5 ms at n = 10³, 5.8 ms at n = 2·10³); extrapolation to n = 10⁴
   gives ≈ 30 ms, far below the 2 s threshold. Optimizing the construction (batching or a hand-written 3×3
   in local coordinates) is not required for performance; the decision at stage 1d
   is driven by the conditioning, not by the time.
2. The functional families (theta, xi1, mu) and the spline evaluation are within the noise
   (ratio 0.8–1.6, the single outlier xi1/H/1000 is 2.4 at an absolute 0.38 ms).
3. The essential result for stage 1d: in the global coordinates of the generating system,
   cond(M_k) on [0,1] grows as n² (≈ 3·10⁴ at n = 10², ≈ 3·10⁶ at n = 10³, ≈ 3·10⁸ at n = 10⁴
   for B; for H it is 4 times larger), so the criterion cond ≤ 10⁸ rejects the basis starting from
   n ≈ 2·10³ (H) and n ≈ 3·10³ (B). In 0.1.0 the det/scale threshold rejected grids with n ≥ 5·10³.
   The move to local coordinates (stage 1d) should remove the dependence of cond on n; until then
   the working range of grid sizes on [0,1] is n ≤ 2·10³ (B) and n ≤ 10³ (H).
