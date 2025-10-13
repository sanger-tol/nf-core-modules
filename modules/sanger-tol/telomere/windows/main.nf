process TELOMERE_WINDOWS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::java-jdk=8.0.112"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/java-jdk:8.0.112--1' :
        'biocontainers/java-jdk:8.0.112--1' }"

    input:
    tuple val(meta), path(telomere)

    output:
    tuple val( meta ), path("*.windows")    , emit: windows
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARNING: This module includes the asmstats script as a module binary in
    // ${moduleDir}/resources/usr/bin/telomere.jar. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or set the option
    // nextflow.enable.moduleBinaries = true in your nextflow.config file.
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def telomere_jar = task.ext.telomere_jar ?: "telomere.jar"
    def telomere_jvm_params = task.ext.telomere_jvm_params ?: "-Xms1g -Xmx1g"
    def telomere_window_cut = task.ext.telomere_window_cut ?: 99.9

    //${projectDir}/bin/

    """
    java ${telomere_jvm_params} \\
        -cp ${telomere_jar} \\
        FindTelomereWindows $telomere $telomere_window_cut \\
        > ${prefix}.windows

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomere: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def telomere = task.ext.telomere ?: ''
    """
    touch ${prefix}.windows

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomere: $VERSION
    END_VERSIONS
    """

}
