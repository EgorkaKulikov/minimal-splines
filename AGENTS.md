# AGENTS.md — правила работы в репозитории `minimal-splines`

Файл обязателен к прочтению перед любым изменением. `CLAUDE.md` ссылается сюда и
второй копии правил не содержит.

## 1. Назначение

Библиотека квадратичных минимальных сплайнов: сетка с тройными краевыми узлами
(`splines.Grid`), порождающие системы `phi = (1, rho, sigma)` (`splines.GeneratingSystem`),
базис `{omega_j}` (`splines.MinimalSplineBasis`), аппроксимационные функционалы и
(квази)проекторы (`splines.functionals.*`), метрика `E_h` (`splines.metrics`).
Публикуется как Maven-артефакт `io.github.egorkakulikov:minimal-splines`.

Библиотека — исследовательский слой: формулы соответствуют публикациям, перечисленным
в `docs/REFERENCES.md`, и каждая конструкция имеет там статус («Подтверждено» /
«Адаптация» / «Без источника»).

## 2. Архитектурные границы

```
numerical-core  <--  minimal-splines  <--  integral-equations (и любые другие потребители)
```

Разрешено импортировать:
- `kotlin.*`, `java.*`;
- `numerics.*` и `numerics.backend.*` из артефакта `numerical-core`
  (`GaussLegendre`, `LinearAlgebra`, `NumericsContext`, `ParallelAssembly`, `orders`, …).

Запрещено:
- `solvers.*`, `problems.*`, `demo.*` — код решателей интегральных уравнений;
- любые понятия уровня интегральных уравнений в коде, KDoc и именах: Фредгольм,
  Вольтерра, Урысон, ядро уравнения, регуляризация, Nyström. Если конструкция нужна
  только решателю (например, стабилизатор Тихонова на сплайнах), её место — в
  `integral-equations`;
- зависимость на `numerical-core` через `project(":...")` или `files("../numerical-core/...")`:
  только опубликованный артефакт. Проверяется задачей `verifyArtifactDependencies`
  (входит в `check`).

Обратная зависимость (`numerical-core -> minimal-splines`) недопустима: если примитив
нужен обеим библиотекам и не знает о сплайнах — он переезжает в `numerical-core`.

## 3. Публичный API и внутренние детали

Публичный API (менять только с обновлением README, KDoc и потребителей):
- `splines`: `Grid` (с фабриками `uniform`, `quasiUniform`, `geometric`, `graded`),
  `GeneratingSystem` (`B`, `H`, `T`), `MinimalSplineBasis`, `ReferenceSplines`,
  `nonDegenerate`, `cross3`, `dot3`, `det3`, `invert3`;
- `splines.functionals`: `ApproxFunctional`, `FunctionalFamily`, `ValueFunctional`,
  `DerivFunctional`, `SecondDerivFunctional`, `ProjFunctionals`, `DeBoorFixFunctionals`,
  `DiscreteDeBoorFixFunctionals`, `AveragingFunctionals`, `ThreePointFunctionals`,
  `SupportPoints`;
- `splines.metrics`: `errorEh`, `DEFAULT_CONTROL_REFINEMENT`.

Внутреннее (`internal`, не расширять видимость):
- `splines/Degeneracy.kt`: `DEGENERACY_RELATIVE_EPS`, `cancellationScale`,
  `isSignificant`, `dot3Scale`, `det3Scale`;
- `MinimalSplineBasis.computeA` (вектор `a_j` аппроксимационного соотношения; нужен
  `AveragingFunctionals`).

Типы `numerical-core` входят в сигнатуры публичного API (`NumericsContext`,
`LinearAlgebra`), поэтому зависимость объявлена `api`. Не делать публичным то, что
раньше было публичным лишь из-за общего модуля.

## 4. Команды

```bash
./gradlew test                      # все тесты, единицы секунд
./gradlew fastTest                  # тег fast (совпадает с test по составу)
./gradlew check                     # test + verifyArtifactDependencies + koverVerify
./gradlew publishToMavenLocal       # артефакт в ~/.m2
../../gradlew -p examples/standalone-consumer run   # контракт стороннего потребителя
```

Бэкенд линейной алгебры в тестах: `multik`; второй прогон —
`./gradlew test -Dnumerics.backend=reference`.

Предварительно `numerical-core` должен быть опубликован в `~/.m2`
(`./gradlew publishToMavenLocal` в его репозитории) либо доступен по
`-PnumericsRepositoryUrl=...`.

## 5. Численные инварианты

Обязаны выполняться на всех четырёх сетках (`uniform`, `quasiUniform`, `graded`,
`geometric`) и всех трёх системах (`B`, `H`, `T`); проверяются
`SplineCoreHealthCheckTest`:
- разбиение единицы: `sum_j omega_j(t) = 1` (следствие `phi_0 == 1`);
- биортогональность проекторов: `theta_i(omega_j) = delta_ij`,
  `xi^<r>_i(omega_j) = delta_ij` для `r = 0, 1, 2`, включая краевые индексы;
- идемпотентность проекторов на сплайнах: `P_chi s = s`;
- точность всех семейств на `span{1, rho, sigma}`;
- совпадение общего построения `M_k^{-1} phi(t)` с замкнутыми формулами
  `ReferenceSplines.omegaB`, `omegaH`; коэффициенты `theta` на равномерной сетке —
  `{1/14, -2/7, 10/7, -2/7, 1/14}`.

Пороги вырожденности ОТНОСИТЕЛЬНЫ масштабу выражения: `DEGENERACY_RELATIVE_EPS = 1e-12`
относительно суммы модулей слагаемых (`cancellationScale`, `det3Scale`, `dot3Scale`).
Абсолютные пороги запрещены — они ломаются на отрезках масштаба `1e-6` и `1e6`
(`DegeneracyScaleTest`). Допуск включения точки разбиения
`Grid.breakpointInclusionEps = 1e-15 * max(1, |b - a|)` относителен только вверх — это
сознательно, ради побитового совпадения с эталонами на `[0, 1]`.

Кратность узла определяется ПО ИНДЕКСУ (`Grid.isCoincident`), не сравнением `Double`
(`GridCoincidenceTest`). Точка вне `[a, b]` — исключение, не экстраполяция.

## 6. Численная нейтральность и эталоны

В этом репозитории эталонных снимков нет: тесты проверяют инварианты. Но сплайны и
функционалы стоят под всеми характеризационными гейтами `integral-equations`
(`baseline-eh.tsv`, `baseline-extra.tsv`, допуск `1e-9`), и любое изменение формул,
порядка суммирования, порогов или путей вычисления здесь ОБЯЗАНО быть проверено там:

```bash
./gradlew publishToMavenLocal
cd ../integral-equations
./gradlew characterizationTest extraCharacterizationTest   # эталоны машинно-зависимы: multik, aarch64
./gradlew fastTest
```

Правила:
- расхождение с эталоном — сначала поиск изменения вычислительного пути, а не правка
  допуска или пересъём эталона;
- пересъём эталона делается только в `integral-equations` по протоколу
  `docs/baseline-changes.md` того репозитория, с указанием причины и коммита библиотеки;
- структурный рефакторинг (перенос, переименование, смена видимости) обязан быть
  нейтрален побитово; изменение алгоритма — отдельный коммит с явной пометкой.

## 7. Regression-тесты

Каждый исправленный дефект получает тест, названный по дефекту, с KDoc: в чём была
ошибка, почему её не ловили существующие тесты, что проверяется теперь. Тест обязан
ПАДАТЬ на коде до исправления (проверить, откатив правку локально). Образцы:
`DegeneracyScaleTest`, `GridCoincidenceTest`, `MinimalSplineBasisExtraTest`.

## 8. Добавление порождающей системы или семейства функционалов

Новая `GeneratingSystem`:
1. первичный источник добавляется в `docs/REFERENCES.md` (раздел 1, таблица и список
   литературы); без источника — статус «Без источника» явно;
2. система включается в `allSystems` теста `SplineCoreHealthCheckTest`: разбиение
   единицы, биортогональность, точность на `span{1, rho, sigma}` на четырёх сетках;
3. при наличии явной формулы — эталон в `ReferenceSplines` и тест совпадения.

Новое семейство функционалов (наследник `FunctionalFamily`):
1. источник в `docs/REFERENCES.md` (раздел 2) — обязателен;
2. корректные `isProjector`, `usesDerivative`, `usesSecondDerivative`;
3. проектор — тест биортогональности и идемпотентности; квазиинтерполянт — тест
   точности на `span{1, rho, sigma}`; краевые `j = -2`, `n - 1` — явно оговорены и
   покрыты;
4. все знаменатели защищены `isSignificant(value, scale)`, а не абсолютным порогом;
5. семейство добавляется в списки `projectorFamilies` / `families` тестов
   `projectorsAreIdempotentOnSplines` и `allFamiliesAreExactOnGeneratingSpan`
   (`SplineCoreHealthCheckTest`);
6. затем — протокол раздела 6 (публикация и гейты `integral-equations`).

Не выдумывать формулы и источники: формула без публикации — «Без источника» в
REFERENCES и явная оговорка в KDoc.

## 9. Язык и стиль

- Комментарии, KDoc, сообщения об ошибках, документация — на русском; идентификаторы —
  на английском; математические обозначения — как в KDoc и REFERENCES (`omega_j`,
  `phi`, `rho`, `sigma`, `theta`, `xi^<r>`, `mu`, `lambda`).
- KDoc объясняет ЗАЧЕМ и ПОЧЕМУ ТАК (обоснование порогов, отвергнутые альтернативы),
  а не пересказывает код.
- Никаких compatibility-обёрток старых пакетов `numerics.functionals` — внешних
  пользователей нет.

## 10. Соседние репозитории

| Репозиторий | Что там | Зависимость |
|---|---|---|
| `numerical-core` | квадратура, линейная алгебра, бэкенды, `NumericsContext`, `orders`/`reliableOrders` | эта библиотека зависит от него (`api`) |
| `integral-equations` | решатели, модельные задачи, демонстрации, характеризационные и верификационные гейты, сверка со SciPy | зависит от этой библиотеки |
| `examples/standalone-consumer` | отдельная Gradle-сборка: сторонний потребитель без `integral-equations` | зависит от обеих библиотек |

Совместное изменение API двух библиотек: сначала `numerical-core` (публикация), затем
здесь (обновить `numericalCoreVersion` в `gradle.properties`), затем `integral-equations`.
