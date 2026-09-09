plugins {
    kotlin("jvm") version "2.0.0"
    `java-library`
    `maven-publish`
    id("org.jetbrains.kotlinx.kover") version "0.8.3"
}

repositories {
    // numerical-core подключается как ОПУБЛИКОВАННЫЙ артефакт, а не как исходники
    // соседнего репозитория. Порядок: mavenLocal первым — локальная сборка
    // (`./gradlew publishToMavenLocal` в numerical-core) имеет приоритет; затем
    // GitHub Packages — оттуда артефакт берёт CI. GitHub Packages требует аутентификацию
    // даже на чтение: переменные окружения GITHUB_ACTOR/GITHUB_TOKEN (в GitHub Actions —
    // встроенный токен) или gpr.user/gpr.token в ~/.gradle/gradle.properties (токен с
    // read:packages). Фильтр content ограничивает этот репозиторий группой библиотеки,
    // чтобы Gradle не ходил в GitHub Packages за остальными зависимостями.
    mavenLocal()
    maven {
        name = "GitHubPackagesNumericalCore"
        url = uri("https://maven.pkg.github.com/EgorkaKulikov/numerical-core")
        credentials {
            username = providers.environmentVariable("GITHUB_ACTOR").orNull ?: providers.gradleProperty("gpr.user").orNull ?: ""
            password = providers.environmentVariable("GITHUB_TOKEN").orNull ?: providers.gradleProperty("gpr.token").orNull ?: ""
        }
        content { includeGroup("io.github.egorkakulikov") }
    }
    mavenCentral()
    // Дополнительный реестр по свойству `numericsRepositoryUrl` (см. README, раздел «Подключение»).
    providers.gradleProperty("numericsRepositoryUrl").orNull?.let { maven(url = uri(it)) }
}

val numericalCoreVersion: String by project

dependencies {
    // `api`, а не `implementation`: типы numerical-core (`GaussLegendre`, `NumericsContext`,
    // `LinearAlgebra`) входят в сигнатуры публичного API этой библиотеки
    // (`FunctionalFamily(basis, ctx)`, `SupportPoints`), поэтому потребитель обязан
    // видеть их на compile classpath транзитивно.
    api("io.github.egorkakulikov:numerical-core:$numericalCoreVersion")

    testImplementation(kotlin("test"))
    testImplementation("org.junit.jupiter:junit-jupiter:5.10.2")
    testRuntimeOnly("org.junit.platform:junit-platform-launcher")
}

kotlin {
    jvmToolchain(21)
}

java {
    withSourcesJar()
}

// Бэкенд линейной алгебры в тестах: см. numerical-core/build.gradle.kts. Здесь он
// нужен семействам функционалов theta/mu/lambda, решающим малые СЛАУ в конструкторе.
val numericsBackend: String = System.getProperty("numerics.backend") ?: "auto"

tasks.test {
    useJUnitPlatform { excludeTags("golden-generate") }
    systemProperty("numerics.backend", numericsBackend)
}

// Формирование эталонов поведения (src/test/resources/golden, см. README.md там же).
// Выполняется только по явному запросу; в `test` и `fastTest` тег golden-generate исключён.
tasks.register<Test>("regenerateGolden") {
    group = "verification"
    description = "Перегенерировать golden-эталоны в src/test/resources/golden"
    testClassesDirs = sourceSets["test"].output.classesDirs
    classpath = sourceSets["test"].runtimeClasspath
    useJUnitPlatform { includeTags("golden-generate") }
    systemProperty("golden.dir", layout.projectDirectory.dir("src/test/resources/golden").asFile.absolutePath)
    systemProperty("golden.version", project.version.toString())
    systemProperty("numerics.backend", numericsBackend)
    outputs.upToDateWhen { false }
}

/** Быстрый набор (тег `fast`); в этой библиотеке совпадает с `test` по составу. */
tasks.register<Test>("fastTest") {
    group = "verification"
    description = "Быстрый набор тестов (тег fast)"
    testClassesDirs = sourceSets["test"].output.classesDirs
    classpath = sourceSets["test"].runtimeClasspath
    useJUnitPlatform { includeTags("fast"); excludeTags("golden-generate") }
    systemProperty("numerics.backend", numericsBackend)
}

kover {
    currentProject {
        instrumentation {
            disabledForTestTasks.add("fastTest")
            disabledForTestTasks.add("regenerateGolden")
        }
    }
}

// --- Проверка независимости от исходников соседних репозиториев -----------------
// Гарантирует, что classpath компиляции содержит numerical-core ТОЛЬКО как jar-артефакт:
// случайная `project(":...")`/`files("../numerical-core/build/...")` зависимость
// провалит задачу. Входит в `check`.
tasks.register("verifyArtifactDependencies") {
    group = "verification"
    description = "Убедиться, что numerical-core подключён как артефакт, а не как исходники"
    val classpath = configurations.compileClasspath
    doLast {
        val offenders = classpath.get().files.filter { f ->
            !f.name.endsWith(".jar") || f.path.contains("${File.separator}numerical-core${File.separator}build${File.separator}")
        }
        check(offenders.isEmpty()) {
            "Зависимости обязаны быть jar-артефактами из репозитория Maven, найдено: $offenders"
        }
        val core = classpath.get().files.filter { it.name.startsWith("numerical-core-") && it.name.endsWith(".jar") }
        check(core.size == 1) { "Ожидался ровно один артефакт numerical-core на classpath, найдено: $core" }
        logger.lifecycle("numerical-core подключён как артефакт: ${core.single().name}")
    }
}
tasks.named("check") { dependsOn("verifyArtifactDependencies") }

// --- Публикация ---------------------------------------------------------------
publishing {
    publications {
        create<MavenPublication>("maven") {
            from(components["java"])
            pom {
                name.set("minimal-splines")
                description.set(
                    "Квадратичные минимальные сплайны на Kotlin/JVM: сетки с кратными узлами, " +
                        "порождающие системы (полиномиальная, гиперболическая, тригонометрическая), " +
                        "базис, аппроксимационные функционалы и квазиинтерполяция.",
                )
                url.set("https://github.com/EgorkaKulikov/minimal-splines")
                licenses {
                    license {
                        name.set("Apache License, Version 2.0")
                        url.set("https://www.apache.org/licenses/LICENSE-2.0")
                    }
                }
                developers {
                    developer {
                        id.set("EgorkaKulikov")
                        name.set("Egor Kulikov")
                    }
                }
                scm {
                    url.set("https://github.com/EgorkaKulikov/minimal-splines")
                }
            }
        }
    }
    repositories {
        // Пример подключения удалённого репозитория (раскомментировать и задать свойства):
        // maven {
        //     name = "remote"
        //     url = uri(providers.gradleProperty("publishUrl").getOrElse(""))
        //     credentials(PasswordCredentials::class) // remoteUsername / remotePassword
        // }
    }
}
