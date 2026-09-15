import kotlinx.kover.gradle.plugin.dsl.CoverageUnit

plugins {
    kotlin("jvm") version "2.0.0"
    `java-library`
    `maven-publish`
    id("org.jetbrains.kotlinx.kover") version "0.8.3"
    id("org.jetbrains.dokka") version "1.9.20"
}

repositories {
    // numerical-core is consumed as a published artifact. Order: mavenLocal first — a local build
    // (`./gradlew publishToMavenLocal` in numerical-core) takes priority; then
    // GitHub Packages — that is where CI takes the artifact from. GitHub Packages requires authentication
    // even for reading: the environment variables GITHUB_ACTOR/GITHUB_TOKEN (in GitHub Actions —
    // the built-in token) or gpr.user/gpr.token in ~/.gradle/gradle.properties (a token with
    // read:packages). The content filter restricts this repository to the group of the library,
    // so that Gradle does not go to GitHub Packages for the remaining dependencies.
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
    // An additional registry given by the `numericsRepositoryUrl` property (see README, section "Usage").
    providers.gradleProperty("numericsRepositoryUrl").orNull?.let { maven(url = uri(it)) }
}

val numericalCoreVersion: String by project

dependencies {
    // `api` rather than `implementation`: numerical-core types (`GaussLegendre`, `NumericsContext`,
    // `LinearAlgebra`) appear in the signatures of the public API of this library
    // (`FunctionalFamily(basis, ctx)`), so a consumer must
    // see them on the compile classpath transitively.
    api("io.github.egorkakulikov:numerical-core:$numericalCoreVersion")

    testImplementation(kotlin("test"))
    testImplementation("net.jqwik:jqwik:1.9.2")
    testImplementation("org.junit.jupiter:junit-jupiter:5.11.4")
    testRuntimeOnly("org.junit.platform:junit-platform-launcher")
}

kotlin {
    jvmToolchain(21)
    explicitApi()
}

java {
    withSourcesJar()
}

// BLAS/LAPACK implementation used in the tests (the numerics.backend property): required by the
// theta/mu/lambda functional families, which solve small linear systems in their constructor.
val numericsBackend: String = System.getProperty("numerics.backend") ?: "auto"

tasks.test {
    useJUnitPlatform { excludeTags("golden-generate") }
    systemProperty("numerics.backend", numericsBackend)
    // Extra stack for the test JVM: the multithreaded dgetrf of the system OpenBLAS on the default
    // thread stack terminated the JVM with a signal (exit 139); see numerical-core.
    jvmArgs("-Xss8m")
}

// Generation of the golden references (src/test/resources/golden, see the README.md there).
// Runs only on explicit request; in `test` and `fastTest` the golden-generate tag is excluded.
tasks.register<Test>("regenerateGolden") {
    group = "verification"
    description = "Regenerate the golden references in src/test/resources/golden"
    testClassesDirs = sourceSets["test"].output.classesDirs
    classpath = sourceSets["test"].runtimeClasspath
    useJUnitPlatform { includeTags("golden-generate") }
    systemProperty("golden.dir", layout.projectDirectory.dir("src/test/resources/golden").asFile.absolutePath)
    systemProperty("golden.version", project.version.toString())
    systemProperty("numerics.backend", numericsBackend)
    outputs.upToDateWhen { false }
}

// Fast suite (the `fast` tag); in this library it has the same content as `test`.
tasks.register<Test>("fastTest") {
    group = "verification"
    description = "Fast test suite (the fast tag)"
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
        sources {
            // The performance measurements are not library code and do not take part in coverage.
            excludedSourceSets.add("benchmark")
        }
    }
    reports {
        filters {
            excludes {
                packages("splines.bench")
            }
        }
        // Coverage bar: enforced by the `koverVerify` task, which is part of `check`.
        // Measured with both BLAS/LAPACK implementations: lines 95.8 %, branches 90.9 %;
        // the threshold is the actual value minus 1 %, so that the result does not depend on the machine.
        verify {
            rule("Line coverage") {
                minBound(94)
            }
            rule("Branch coverage") {
                bound {
                    minValue = 89
                    coverageUnits = CoverageUnit.BRANCH
                }
            }
        }
    }
}

tasks.check {
    dependsOn("koverVerify")
}

// --- Documentation ------------------------------------------------------------
// HTML documentation of the public API: ./gradlew dokkaHtml (the result is in build/dokka/html).
// Public symbols without KDoc are reported as warnings.
tasks.dokkaHtml {
    moduleName.set("minimal-splines")
    dokkaSourceSets.configureEach {
        includeNonPublic.set(false)
        reportUndocumented.set(true)
        jdkVersion.set(21)
    }
}

// --- Publication --------------------------------------------------------------
publishing {
    publications {
        create<MavenPublication>("maven") {
            from(components["java"])
            pom {
                name.set("minimal-splines")
                description.set(
                    "A library of quadratic minimal splines and quasi-interpolation for the JVM platform: " +
                        "a basis built from an arbitrary generating system (polynomial, hyperbolic, " +
                        "trigonometric) and local approximation functionals without solving a global linear system.",
                )
                url.set("https://github.com/EgorkaKulikov/minimal-splines")
                inceptionYear.set("2026")
                licenses {
                    license {
                        name.set("Apache License, Version 2.0")
                        url.set("https://www.apache.org/licenses/LICENSE-2.0")
                        distribution.set("repo")
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
                    connection.set("scm:git:https://github.com/EgorkaKulikov/minimal-splines.git")
                    developerConnection.set("scm:git:ssh://git@github.com/EgorkaKulikov/minimal-splines.git")
                }
            }
        }
    }
    repositories {
        // Publication to GitHub Packages is performed from CI on a version tag by the task
        // publishAllPublicationsToGitHubPackagesRepository; the credentials are GITHUB_ACTOR/GITHUB_TOKEN
        // or gpr.user/gpr.token. Without them the repository is declared but unreachable; publishToMavenLocal
        // does not depend on it.
        maven {
            name = "GitHubPackages"
            url = uri("https://maven.pkg.github.com/EgorkaKulikov/minimal-splines")
            credentials {
                username = providers.environmentVariable("GITHUB_ACTOR").orNull ?: providers.gradleProperty("gpr.user").orNull
                password = providers.environmentVariable("GITHUB_TOKEN").orNull ?: providers.gradleProperty("gpr.token").orNull
            }
        }
    }
}

// Performance measurements of the public API; they are part of neither the artifact nor test.
val benchmark: SourceSet by sourceSets.creating {
    compileClasspath += sourceSets["main"].output + configurations["runtimeClasspath"]
    runtimeClasspath += output + compileClasspath
}
tasks.register<JavaExec>("benchmark") {
    description = "Performance measurements of basis and functional construction; grid sizes are given via -Pbench.args=\"100 1000 10000\""
    group = "verification"
    classpath = benchmark.runtimeClasspath
    mainClass.set("splines.bench.BenchKt")
    maxHeapSize = "4g"
    args = (project.findProperty("bench.args") as String? ?: "100 1000 10000").split(" ").filter { it.isNotBlank() }
    systemProperty("numerics.backend", numericsBackend)
}
