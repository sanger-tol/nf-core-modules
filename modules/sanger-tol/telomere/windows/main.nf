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
    // WARNING: This module includes the telomere.jar binary as a module binary in
    // ${moduleDir}/resources/usr/bin/telomere.jar. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or set the option
    // nextflow.enable.moduleBinaries = true in your nextflow.config file.

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "TELOMERE_WINDOWS module does not support Conda. Please use Docker / Singularity instead."
    }

    def prefix      = task.ext.prefix ?: "${meta.id}"
    def args        = task.ext.args   ?: ""
    def VERSION     = "1.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    // Dynamically generate java mem needs based on task.memory
    // Taken from: nf-core/umicollapse
    def max_heap_size_mega = (task.memory.toMega() * 0.9).intValue()
    def max_stack_size_mega = 999 //most java jdks will not allow Xss > 1GB, so fixing this to the allowed max

    // Use groovy to move through list of expected jar locations
    // Error if not found.
    def jar_locations = ["${moduleDir}/resources/usr/bin/telomere.jar", "${projectDir}/bin/telomere.jar", task.ext.jar]
    def jar = jar_locations.find { file(it).exists() }
    if(!jar) {
        log.error("ERROR: Could not locate a telomere JAR file!")
    }

    """
    java \\
        -Xmx${max_heap_size_mega}M \\
        -Xss${max_stack_size_mega}M \\
        -cp ${jar} \\
        FindTelomereWindows $telomere \\
        $args \\
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
