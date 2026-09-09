# minimal-splines

Библиотека квадратичных минимальных сплайнов и квазиинтерполяции на Kotlin/JVM:
сетка с кратными узлами, порождающие системы, базис `{omega_j}`, аппроксимационные
функционалы `theta`, `xi`, `xitilde`, `mu`, `lambda` и метрика ошибки `E_h`.

Библиотека — второй слой в цепочке репозиториев:

```
numerical-core  <--  minimal-splines  <--  integral-equations
```

`minimal-splines` зависит только от [`numerical-core`](https://github.com/EgorkaKulikov/numerical-core)
(квадратура Гаусса–Лежандра, линейная алгебра, `NumericsContext`) и ничего не знает об
интегральных уравнениях. Любой проект, которому нужны минимальные сплайны, подключает
эту библиотеку напрямую — без `integral-equations`; исполняемое доказательство —
[`examples/standalone-consumer`](examples/standalone-consumer/README.md).

## Что входит

| Объект | Класс | Назначение |
|---|---|---|
| Сетка `X` | `splines.Grid` | `a = x_{-2} = x_{-1} = x_0 < x_1 < … < x_n = x_{n+1} = x_{n+2} = b`: узлы кратности 3 на концах, `x(j)` для `j = -2..n+2`, шаг `h = max (x_{j+1} - x_j)`; фабрики `uniform`, `quasiUniform`, `geometric`, `graded` |
| Порождающая система `phi` | `splines.GeneratingSystem` | вектор-функция `phi(t) = (1, rho(t), sigma(t))` с производными до второго порядка; готовые `B = (1, t, t^2)`, `H = (1, sinh t, cosh t)`, `T = (1, sin t, cos t)` |
| Базис | `splines.MinimalSplineBasis` | сплайны `omega_j`, `j = -2..n-1`, с носителями `[x_j, x_{j+3}]`: `omega`, `omegaDeriv`, `omegaDeriv2`, `activeOmega`, `evalSpline`, `evalSplineDeriv`, `evalSplineDeriv2`, `interval` |
| Замкнутые формулы | `splines.ReferenceSplines` | `omegaB`, `omegaBDeriv`, `omegaH` — независимые от общего построения эталоны для тестов |
| Функционалы | `splines.functionals.*` | `ApproxFunctional`, `FunctionalFamily` и пять семейств (таблица ниже) |
| Опорные точки | `splines.functionals.SupportPoints` | объединение точек семейства `ValueFunctional` с явной индексацией `(j, q) -> r` |
| Метрика | `splines.metrics.errorEh` | `E_h = max |u*(t) - u_h(t)|` на контрольной сетке `refinement * n + 1` точек |

## Что намеренно НЕ входит

- интегральные операторы и решатели (Фредгольм, Вольтерра, Урысон), модельные задачи,
  таблицы сходимости решателей, регуляризация — репозиторий `integral-equations`;
- квадратура, плотная линейная алгебра, бэкенды, обусловленность, `NumericsContext`,
  порядки сходимости `orders`/`reliableOrders` — репозиторий `numerical-core`
  (подключается транзитивно);
- сплайны степени, отличной от 2, и сетки без тройных краевых узлов.

## Подключение

Координаты: `io.github.egorkakulikov:minimal-splines:0.1.0`. Зависимость на
`numerical-core` объявлена как `api` и приходит транзитивно.

```kotlin
repositories {
    mavenCentral()
    mavenLocal() // после `publishToMavenLocal` в обеих библиотеках
    // numerical-core также доступен из GitHub Packages — см. ниже
}

dependencies {
    implementation("io.github.egorkakulikov:minimal-splines:0.1.0")
}
```

Локальная публикация обеих библиотек (по порядку):

```bash
(cd ../numerical-core && ./gradlew publishToMavenLocal)
./gradlew publishToMavenLocal
```

`numerical-core` публикуется в GitHub Packages
(`https://maven.pkg.github.com/EgorkaKulikov/numerical-core`) — так его берёт CI, так он
подключён в `build.gradle.kts` (после `mavenLocal()`). Для чтения нужен токен с правом
`read:packages`: переменные окружения `GITHUB_ACTOR`/`GITHUB_TOKEN` или `gpr.user`/`gpr.token`
в `~/.gradle/gradle.properties`; локально можно обойтись `publishToMavenLocal`. Сам
`minimal-splines` в удалённый реестр пока не публикуется: точка настройки — блок
`publishing.repositories` в `build.gradle.kts` и свойство `numericsRepositoryUrl` у потребителей.

## Минимальный пример

Проекция функции на пространство сплайнов и оценка погрешности (тот же сценарий, что в
`MinimalSplineBasisExtraTest.evalSplineReproducesLinear`):

```kotlin
import splines.GeneratingSystem
import splines.Grid
import splines.MinimalSplineBasis
import splines.functionals.ProjFunctionals
import splines.metrics.errorEh

fun main() {
    val grid = Grid.uniform(n = 16, a = 0.0, b = 1.0)
    val basis = MinimalSplineBasis(GeneratingSystem.B, grid)
    val theta = ProjFunctionals(basis)

    val f = { t: Double -> Math.exp(Math.sin(3.0 * t)) }
    val c = theta.projectorCoeffs(f)                       // chi_j(f), j = -2..n-1
    val spline = { t: Double -> basis.evalSpline(c, t) }   // P_theta f = sum_j chi_j(f) omega_j

    println("E_h = ${errorEh(f, spline, grid)}")            // ~2e-4 при n = 16, порядок 3
    println("f'(0.5) ~ ${basis.evalSplineDeriv(c, 0.5)}")
}
```

Для семейства `xi` (де Бура–Фикса) в `projectorCoeffs` передаются также производные
`gD` и, при `r = 0`, `gDD`.

## Основные понятия API

- **`Grid`.** Индексация математическая: `x(j)` для `j = -2..n+2`, `idx(j) = j + 2`.
  Внутренние узлы строго возрастают (проверяется в конструкторе); кратность краевых
  узлов задаётся паддингом, поэтому `isCoincident(j)` отвечает по индексу, а не сравнением
  `Double`. Массив `breakpoints` (`x_0..x_n`) — горячие данные, read-only по соглашению.
- **`GeneratingSystem`.** Требуется `rho`, `sigma` и их производные до второго порядка;
  `wronskian(t) = det(phi, phi', phi'')` — проверка невырожденности.
- **`MinimalSplineBasis`.** На `(x_k, x_{k+1})` активны три сплайна
  `omega_{k-2}, omega_{k-1}, omega_k`; значения получаются как `M_k^{-1} phi(t)`, где
  `M_k = (a_{k-2} | a_{k-1} | a_k)` — матрица аппроксимационного соотношения.
  Точка вне `[a, b]` — ошибка, а не экстраполяция.
- **`FunctionalFamily`.** Семейство `{chi_j}_{j=-2}^{n-1}` и (квази)проектор
  `P_chi g = sum_j chi_j(g) omega_j`. Флаги `isProjector` (`P^2 = P`, биортогональность),
  `usesDerivative`, `usesSecondDerivative` сообщают, какие производные `g` нужны.
  `projectorCoeffs(g, gD, gDD)` возвращает вектор коэффициентов размера `n + 2`;
  `cChi()` — максимум `absSum()` по семейству (оценка усиления возмущения данных —
  только для семейств из функционалов-значений).
- **`NumericsContext`** (из `numerical-core`). Семейства `theta`, `mu`, `lambda`
  решают малые СЛАУ в конструкторе и принимают контекст последним параметром; по
  умолчанию — `NumericsContext.default()`.

## Семейства функционалов

| Семейство | Класс | Проектор | Нужны производные `g` | Источник |
|---|---|---|---|---|
| `theta` — проекционные, 5-точечная закрытая форма | `ProjFunctionals` | да | нет | [Kulikov, Makarov 2025] |
| `xi^<r>`, `r ∈ {0, 1, 2}` — де Бура–Фикса | `DeBoorFixFunctionals(basis, r)` | да | `g'`; при `r = 0` также `g''` | [Kulikov, Makarov 2019b] |
| `xitilde^<r>`, `r ∈ {1, 2}` — дискретизованные де Бура–Фикса | `DiscreteDeBoorFixFunctionals(basis, r)` | нет | нет | собственная конструкция |
| `mu` — усредняющие по вспомогательной сетке `y_j = x_{j+1} + theta (x_{j+2} - x_{j+1})` | `AveragingFunctionals(basis, theta = 0.5)` | нет | нет | [Kulikov, Makarov 2022] |
| `lambda` — трёхточечные по `x_{j+1}, x_{j+3/2}, x_{j+2}` | `ThreePointFunctionals(basis, thetaHat = 0.5)` | нет | нет | [Kulikov, Makarov 2022] |

Пять классов дают восемь именованных вариантов. Ограничения:

- `xitilde^<0>` не реализован по построению: он опирался бы на `f''` в левом конце
  носителя, где сплайн и его первая производная обращаются в ноль, и замена второй
  производной разностью не воспроизводит порождающую систему даже при `h -> 0`.
- Семейства `xi` требуют производных образа; `xitilde` заменяют производную центральной
  разделённой разностью и работают только со значениями, но перестают быть проекторами.
- Краевые функционалы `xi` при `j = -2`, `n - 1` взяты как чистые значения `u(x_0)`,
  `u(x_n)` — допущение реализации, закрытое тестом биортогональности
  (см. `docs/REFERENCES.md`, раздел 2).
- `ProjFunctionals` строит `theta_j` локальной биортогонализацией (решение СЛАУ
  `theta_j(omega_i) = delta_ij` по узлам и серединам интервалов) на любой сетке; явная
  пятиточечная формула `closedFormInternal` служит независимой сверкой и совпадает с
  построением на равномерной сетке (`SplineCoreHealthCheckTest.closedFormIsUniformGridOnly`).
- Сетки: `uniform`, `quasiUniform` (`Psi(u) = u + amp sin 2 pi u`), `geometric`
  (отношение крайних шагов `R`), `graded` (чередование шагов `s, ratio * s`). Минимальные
  сплайны существуют на локально квазиравномерных сетках; строгая монотонность узлов
  проверяется конструктором `Grid`.

## Численная верификация

Тесты проверяют математические следствия формул, а не совпадение с эталонными числами;
характеризационных эталонов в этой библиотеке нет.

| Свойство | Тест |
|---|---|
| Совпадение общего построения `M_k^{-1} phi(t)` с замкнутыми формулами `omegaB`, `omegaH` | `SplineCoreHealthCheckTest.polynomialSplineMatchesClosedFormB`, `hyperbolicSplineMatchesClosedFormH` |
| Разбиение единицы `sum_j omega_j(t) = 1` | `SplineCoreHealthCheckTest.basisFormsPartitionOfUnity` |
| Биортогональность `chi_i(omega_j) = delta_ij` для `theta` и всех `xi^<r>`, включая краевые индексы | `projectionFunctionalsAreBiorthogonal`, `deBoorFixFunctionalsAreBiorthogonalForAllOrders` |
| Идемпотентность проекторов на сплайнах | `projectorsAreIdempotentOnSplines` |
| Точность всех семейств на `span{1, rho, sigma}` для систем `B`, `H`, `T` | `allFamiliesAreExactOnGeneratingSpan` |
| Опубликованные коэффициенты `theta` `{1/14, -2/7, 10/7, -2/7, 1/14}` | `closedFormCoefficientsMatchPublishedValues` |
| Сетки `uniform`, `quasiUniform`, `graded`, `geometric` во всех проверках выше | `SplineCoreHealthCheckTest` (четыре сетки на каждую проверку) |
| Индексный критерий кратного узла побитово совпадает со сравнением значений | `GridCoincidenceTest` |
| Масштабная инвариантность порогов вырожденности (`[0, 1e-6]`, `[0, 1e6]`) | `DegeneracyScaleTest` |
| `omegaDeriv` согласуется с численной производной; вырожденность на правом тройном крае | `MinimalSplineBasisExtraTest` |
| Объединение опорных точек по допуску, порядок `byAscendingValue`/`byFirstOccurrence` | `SupportPointsTest` |

Внешняя сверка базиса `B` со `scipy.interpolate.BSpline` и сквозные характеризационные
гейты (значения `E_h` решателей против эталонов с допуском `1e-9`) выполняются в
репозитории `integral-equations` против опубликованного артефакта этой библиотеки.
Поэтому любое изменение формул здесь завершается прогоном
`characterizationTest extraCharacterizationTest` там (см. `AGENTS.md`).

## Тесты и сборка

```bash
./gradlew test        # все тесты (единицы секунд)
./gradlew fastTest    # тег fast — совпадает с test по составу
./gradlew check       # test + verifyArtifactDependencies + koverVerify
./gradlew build
```

`verifyArtifactDependencies` проверяет, что `numerical-core` присутствует на classpath
только как jar из репозитория Maven, а не как исходники соседнего каталога.

Бэкенд линейной алгебры в тестах фиксирован: `multik`; переопределение —
`-Dnumerics.backend=reference`.

## Публикация

```bash
./gradlew publishToMavenLocal              # ~/.m2/repository/io/github/egorkakulikov/minimal-splines/0.1.0/
./gradlew publishToMavenLocal -Pversion=0.2.0
```

Версия и группа — в `gradle.properties`; версия `numerical-core`, против которой
собирается библиотека, — свойство `numericalCoreVersion` там же.

## Теоретические источники

Соответствие «элемент реализации → публикация → статус сверки» — в
[`docs/REFERENCES.md`](docs/REFERENCES.md): аппроксимационные соотношения и порождающие
системы (раздел 1), функционалы (раздел 2). Конструкции без опубликованного источника
(`xitilde`, краевые `xi`) помечены там явно.

## Соседние репозитории

| Репозиторий | Роль | Направление зависимости |
|---|---|---|
| `numerical-core` | квадратура, линейная алгебра, бэкенды, контекст, порядки сходимости | `minimal-splines -> numerical-core` |
| `integral-equations` | решатели Фредгольма, Вольтерры, Урысона; характеризационные и верификационные гейты | `integral-equations -> minimal-splines` |
| `examples/standalone-consumer` | образец независимого потребителя без `integral-equations` | `consumer -> {numerical-core, minimal-splines}` |

Обратных зависимостей нет: код этой библиотеки не импортирует `solvers.*`, `problems.*`.

## Участие в разработке

Правила для людей и агентов — в [`AGENTS.md`](AGENTS.md): границы импорта, публичный
API, численные инварианты, обязательные тесты для новых семейств функционалов и
порождающих систем, протокол проверки численной нейтральности через `integral-equations`.

## Лицензия

Apache License 2.0 — см. [`LICENSE`](LICENSE).
