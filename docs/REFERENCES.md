# Источники методов

В документе зафиксировано происхождение каждого метода библиотеки: математической конструкции,
формулы или инженерного решения. Ссылки даны на издания, по которым сверялись формулы; для
элементов без внешнего источника указано место, где приведено обоснование.

## Условные обозначения статуса

- **Источник** — конструкция описана в указанной работе; реализация сверена с формулами
  источника и тестами.
- **Собственная реализация** — инженерное решение библиотеки без внешнего источника;
  обоснование приведено в KDoc указанного элемента и в `docs/ТОЧНОСТЬ.md`.

## 1. Порождающие системы и базис

| Элемент | Реализация | Источник | Статус |
|---|---|---|---|
| Порождающая система φ = (1, ρ, σ) с невырожденным вронскианом; полиномиальная система (1, t, t²) | `GeneratingSystem`, `GeneratingSystem.B` | [Демьянович 1994]; [Makarov 2012] | Источник |
| Гиперболическая система (1, sinh t, cosh t) | `GeneratingSystem.H` | [Kulikov, Makarov 2019a] | Источник |
| Тригонометрическая система (1, sin t, cos t) | `GeneratingSystem.T` | [Kosogorov, Makarov 2017] | Источник |
| Сетка с тройными узлами на концах отрезка | `Grid` | [Kulikov, Makarov 2020] | Источник |
| Аппроксимационное соотношение Σ_j a_j ω_j = φ, векторы a_j по значениям φ, φ' в узлах, базис ω(t) = M_k⁻¹ φ(t) | `MinimalSplineBasis`, `MinimalSplineBasis.computeA` | [Демьянович 1994]; [Makarov 2012] | Источник |
| Явные формулы базисных функций для систем B и H (оракул в тестах) | `src/test/kotlin/splines/ReferenceSplines.kt` | [Makarov 2012]; [Kulikov, Makarov 2019a] | Источник |
| Совпадение с квадратичными B-сплайнами для системы B; рекурсия Кокса–де Бура в тестах | `src/test/kotlin/splines/DeBoorOracleTest.kt` | [de Boor 2001], гл. IX | Источник |
| Локальные координаты интервала ψ_k = T_k φ; ковариантность a_j^{Tφ} = T a_j^{φ} | `LocalFrame`, `GeneratingSystem.localFrame`, `MinimalSplineBasis.localApproximationMatrix` | обоснование в `docs/ТОЧНОСТЬ.md`, раздел «Локальные координаты интервала» | Собственная реализация |
| Критерий вырожденности по числу обусловленности T_k M_k и относительный порог значимости знаменателей | `MinimalSplineBasis.MAX_CONDITION`, `DEGENERACY_RELATIVE_EPS`, `nonDegenerate` | обоснование в `docs/ТОЧНОСТЬ.md` | Собственная реализация |

## 2. Аппроксимационные функционалы

| Элемент | Реализация | Источник | Статус |
|---|---|---|---|
| Проекционные функционалы θ_j: пятиточечная замкнутая форма по узлам x_j, x_{j+3} и серединам интервалов носителя; биортогональность θ_i(ω_j) = δ_ij | `ProjFunctionals` | [Kulikov, Makarov 2025] | Источник |
| Функционалы де Бура–Фикса ξ^⟨r⟩, r ∈ {0, 1, 2}, для полиномиальных сплайнов | `DeBoorFixFunctionals` | [de Boor, Fix 1973] | Источник |
| Функционалы ξ^⟨r⟩ для минимальных сплайнов по произвольной порождающей системе; биортогональность ξ_i(ω_j) = δ_ij | `DeBoorFixFunctionals` | [Kulikov, Makarov 2019b] | Источник |
| Краевые функционалы θ, ξ, ξ̃ при j = −2 и j = n − 1 как значения g(x₀), g(x_n) | `ProjFunctionals`, `DeBoorFixFunctionals`, `DiscreteDeBoorFixFunctionals` | биортогональность проверяется тестами для всех r и всех систем | Собственная реализация |
| Дискретные функционалы де Бура–Фикса ξ̃^⟨r⟩, r ∈ {1, 2}: замена производной разделённой разностью по соседним узлам | `DiscreteDeBoorFixFunctionals` | свойства в `docs/ТОЧНОСТЬ.md`, раздел «Функционалы ξ̃» | Собственная реализация |
| Усредняющие функционалы μ_j по значениям в точках вспомогательной сетки y_j = x_{j+1} + θ (x_{j+2} − x_{j+1}) | `AveragingFunctionals` | [Kulikov, Makarov 2022] | Источник |
| Трёхточечные функционалы λ_j по значениям в x_{j+1}, x_{j+3/2}, x_{j+2} | `ThreePointFunctionals` | [Kulikov, Makarov 2022] | Источник |
| Вычисление коэффициентов функционалов через локальное представление ψ_k вместо φ | `FunctionalFamily` и наследники | обоснование в `docs/ТОЧНОСТЬ.md` | Собственная реализация |

## 3. Метрики

| Элемент | Реализация | Источник | Статус |
|---|---|---|---|
| Погрешность E_h как максимум |g − P_χ g| по равномерной контрольной сетке с заданным числом точек на интервал | `errorEh`, `DEFAULT_CONTROL_REFINEMENT` | стандартное определение дискретной равномерной нормы | Источник |
| Наблюдаемый порядок сходимости log₂(E_h / E_{h/2}) | тесты сходимости; функции `orders`, `reliableOrders` библиотеки numerical-core | стандартное определение | Источник |

## Список литературы

1. **[Демьянович 1994]** Демьянович Ю. К. Локальная аппроксимация на многообразии и минимальные
   сплайны. — СПб.: Изд-во С.-Петерб. ун-та, 1994. — 356 с.
2. **[de Boor, Fix 1973]** de Boor C., Fix G. J. Spline approximation by quasiinterpolants //
   Journal of Approximation Theory. — 1973. — Vol. 8, No. 1. — P. 19–45.
3. **[de Boor 2001]** de Boor C. A Practical Guide to Splines. — Revised ed. — New York:
   Springer, 2001. — (Applied Mathematical Sciences; vol. 27).
4. **[Makarov 2012]** Makarov A. A. Construction of Splines of Maximal Smoothness // Journal of
   Mathematical Sciences. — 2012. — Vol. 178, No. 6. — P. 589–604.
5. **[Kosogorov, Makarov 2017]** Kosogorov O., Makarov A. On Some Piecewise Quadratic Spline
   Functions // Numerical Analysis and Its Applications. — Cham: Springer, 2017. — (Lecture Notes
   in Computer Science; vol. 10187). — P. 448–455.
6. **[Kulikov, Makarov 2019a]** Kulikov E. K., Makarov A. A. On Approximation by Hyperbolic
   Splines // Journal of Mathematical Sciences. — 2019. — Vol. 240, No. 6. — P. 822–832.
7. **[Kulikov, Makarov 2019b]** Kulikov E. K., Makarov A. A. On de Boor–Fix Type Functionals for
   Minimal Splines // Topics in Classical and Modern Analysis. — Cham: Springer, 2019. —
   (Applied and Numerical Harmonic Analysis). — P. 211–225.
8. **[Kulikov, Makarov 2020]** Kulikov E. K., Makarov A. A. Quadratic Minimal Splines with
   Multiple Nodes // Journal of Mathematical Sciences. — 2020. — Vol. 249, No. 2. — P. 256–262.
9. **[Kulikov, Makarov 2022]** Kulikov E. K., Makarov A. A. Construction of Approximation
   Functionals for Minimal Splines // Journal of Mathematical Sciences. — 2022. — Vol. 262,
   No. 1. — P. 84–98.
10. **[Kulikov, Makarov 2025]** Kulikov E. K., Makarov A. A. On Projection-Type Approximation
    Functionals for Minimal Splines // Записки научных семинаров ПОМИ. — 2025. — Т. 542. —
    С. 126–143.
