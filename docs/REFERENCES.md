# Источники: минимальные сплайны и аппроксимационные функционалы

Документ относится к библиотеке `minimal-splines`. Источники схем решения интегральных
уравнений (Фредгольм, Вольтерра, Урысон), регуляризации и независимой верификации
решателей ведутся в репозитории `integral-equations` (`docs/REFERENCES.md`), который
ссылается на настоящий документ за разделами 1–2.

Документ связывает каждый реализованный численный метод с первичной публикацией и
фиксирует статус соответствия: подтверждено сверкой формул, адаптация или
собственная конструкция без опубликованного аналога.

Проверка проводилась сопоставлением формул кода с формулами источников. Для методов,
у которых опубликованного первоисточника нет, это указано явно — такие места **не
следует** сопровождать ссылками на литературу.

## Условные обозначения статуса

| Статус | Значение |
|---|---|
| **Подтверждено** | формула кода дословно совпадает с формулой публикации |
| **Адаптация** | метод перенесён на класс задач, для которого он в источнике не доказан |
| **Без источника** | опубликованного аналога не найдено; конструкция собственная |

---

## 1. Минимальные сплайны и порождающие системы

| Элемент | Реализация | Источник | Статус |
|---|---|---|---|
| Аппроксимационные соотношения, построение `a_j`, `M_k^{-1} phi(t)` | `splines/MinimalSplineBasis.kt` | [Makarov 2012], [Демьянович 1994] | Подтверждено |
| Полиномиальная система `phi^B = (1, t, t^2)` | `splines/GeneratingSystem.kt` | [Makarov 2012], [Демьянович 1994] | Подтверждено |
| Гиперболическая система `phi^H = (1, sinh t, cosh t)` | `splines/GeneratingSystem.kt` | [Kulikov, Makarov 2019a] | Подтверждено |
| Тригонометрическая система `phi^T = (1, sin t, cos t)` | `splines/GeneratingSystem.kt` | [Kosogorov, Makarov 2017] | Подтверждено |
| Явные формулы `ReferenceSplines.omegaB` | `splines/MinimalSplineBasis.kt` | [Makarov 2012] | Подтверждено |
| Явные формулы `ReferenceSplines.omegaH` | `splines/MinimalSplineBasis.kt` | [Kulikov, Makarov 2019a] | Подтверждено (выведено из общей формулы) |
| Узлы кратности 3 на концах отрезка | `splines/Grid.kt` | [Kulikov, Makarov 2020] | Подтверждено |

## 2. Аппроксимационные функционалы

| Элемент | Реализация | Источник | Статус |
|---|---|---|---|
| Проекционные `theta_j` (5-точечная закрытая форма) | `splines/functionals/Functionals.kt`, `ProjFunctionals` | [Kulikov, Makarov 2025] | Подтверждено |
| Функционалы де Бура–Фикса `xi^<0>`, `xi^<1>`, `xi^<2>` | `DeBoorFixFunctionals` | [Kulikov, Makarov 2019b] | Подтверждено (послагаемо) |
| Биортогональность `xi_i(omega_j) = delta_ij` | `DeBoorFixFunctionals` | [Kulikov, Makarov 2019b] | Подтверждено |
| Краевые `xi` при `j = -2`, `n-1` (чистые значения) | `DeBoorFixFunctionals.buildXi` | — | **Без источника**: допущение реализации, закрытое тестом биортогональности `SplineCoreHealthCheckTest.deBoorFixFunctionalsAreBiorthogonalForAllOrders` (включает краевые индексы для всех `r` и всех базисов) |
| Усредняющие `mu_j` | `AveragingFunctionals` | [Kulikov, Makarov 2022] | Подтверждено |
| Трёхточечные `lambda_j` | `ThreePointFunctionals` | [Kulikov, Makarov 2022] | Подтверждено |
| Дискретизованные `xitilde^<1>`, `xitilde^<2>` | `DiscreteDeBoorFixFunctionals` | — | **Без источника**: новизна работы, в литературе ранее не рассматривались |

## 3. Методы, использующие сплайны и функционалы

Схемы решения интегральных уравнений второго и первого рода (коллокация, итерация
Слоана, Кулкарни, Nyström, комбинированный Nyström, регуляризация), нелинейное
уравнение Урысона и таблица «метод → средства верификации» описаны в репозитории
`integral-equations`, файл `docs/REFERENCES.md`, разделы 3–6. Там же зафиксированы
статусы «Адаптация» и «Без источника» для конструкций уровня решателей.

Из средств независимой верификации к настоящей библиотеке относятся два:

| Средство | Реализация здесь | Степень независимости |
|---|---|---|
| Математические инварианты (биортогональность, разбиение единицы, идемпотентность, точность на `span{1, rho, sigma}`) | `src/test/kotlin/splines/SplineCoreHealthCheckTest.kt` | Высокая: свойства выведены из теории, а не из реализации |
| Сверка с явными замкнутыми формулами (`ReferenceSplines.omegaB`, `omegaH`; `ProjFunctionals.closedFormInternal`) | `SplineCoreHealthCheckTest.kt` | Высокая: формулы выписаны независимо от общего построения `M_k^{-1} phi(t)` |

Внешняя сверка базиса `B` со `scipy.interpolate.BSpline` выполняется в
`integral-equations` (задача `scipyVerify`, скрипт `tools/verify_with_scipy.py`), потому
что использует общий с решателями механизм выгрузки артефактов.

---

## Список литературы

1. **[Демьянович 1994]** Демьянович Ю. К. *Локальная аппроксимация на многообразии и
   минимальные сплайны.* — СПб.: Изд-во С.-Петерб. ун-та, 1994. — 356 с.

2. **[Makarov 2012]** Makarov A. A. Construction of Splines of Maximal Smoothness //
   Journal of Mathematical Sciences. — 2012. — Vol. 178, No. 6. — P. 589–604.

3. **[Kosogorov, Makarov 2017]** Kosogorov O., Makarov A. On Some Piecewise Quadratic
   Spline Functions // Numerical Analysis and Its Applications. Lecture Notes in
   Computer Science, Vol. 10187. — 2017. — P. 448–455.

4. **[Kulikov, Makarov 2019a]** Kulikov E. K., Makarov A. A. On Approximation by
   Hyperbolic Splines // Journal of Mathematical Sciences. — 2019. — Vol. 240, No. 6. —
   P. 822–832.

5. **[Kulikov, Makarov 2019b]** Kulikov E. K., Makarov A. A. On de Boor–Fix Type
   Functionals for Minimal Splines // Topics in Classical and Modern Analysis (Applied
   and Numerical Harmonic Analysis). — Springer, 2019. — P. 211–225.

6. **[Kulikov, Makarov 2020]** Kulikov E. K., Makarov A. A. Quadratic minimal splines
   with multiple nodes // Journal of Mathematical Sciences. — 2020. — Vol. 249, No. 2. —
   P. 256–262.

7. **[Kulikov, Makarov 2022]** Kulikov E. K., Makarov A. A. Construction of
   Approximation Functionals for Minimal Splines // Journal of Mathematical Sciences. —
   2022. — Vol. 262, No. 1. — P. 84–98.

8. **[Kulikov, Makarov 2025]** Kulikov E. K., Makarov A. A. On Projection-Type
    Approximation Functionals for Minimal Splines // Записки научных семинаров ПОМИ. —
    2025. — Т. 542. — С. 126–143.
