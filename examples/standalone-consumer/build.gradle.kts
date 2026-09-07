plugins {
    kotlin("jvm") version "2.0.0"
    application
}

repositories {
    mavenCentral()
    // Артефакты библиотек: локальный репозиторий Maven (после `publishToMavenLocal`
    // в numerical-core и minimal-splines) либо удалённый из `-PnumericsRepositoryUrl=...`.
    mavenLocal()
    providers.gradleProperty("numericsRepositoryUrl").orNull?.let { maven(url = uri(it)) }
}

dependencies {
    implementation("io.github.egorkakulikov:numerical-core:0.1.0")
    implementation("io.github.egorkakulikov:minimal-splines:0.1.0")
    // integral-equations здесь НЕТ и быть не должно: это контракт независимости библиотек.
}

kotlin {
    jvmToolchain(21)
}

application {
    mainClass.set("consumer.MainKt")
}
