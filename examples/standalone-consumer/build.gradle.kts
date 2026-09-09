plugins {
    kotlin("jvm") version "2.0.0"
    application
}

repositories {
    // Артефакты библиотек: локальный репозиторий Maven (после `publishToMavenLocal`
    // в numerical-core и minimal-splines) — первым; numerical-core дополнительно доступен
    // из GitHub Packages (нужны GITHUB_ACTOR/GITHUB_TOKEN или gpr.user/gpr.token);
    // остальное — из mavenCentral либо реестра из `-PnumericsRepositoryUrl=...`.
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
    providers.gradleProperty("numericsRepositoryUrl").orNull?.let { maven(url = uri(it)) }
}

dependencies {
    implementation("io.github.egorkakulikov:numerical-core:1.0.0")
    implementation("io.github.egorkakulikov:minimal-splines:0.1.0")
    // integral-equations здесь НЕТ и быть не должно: это контракт независимости библиотек.
}

kotlin {
    jvmToolchain(21)
}

application {
    mainClass.set("consumer.MainKt")
}
