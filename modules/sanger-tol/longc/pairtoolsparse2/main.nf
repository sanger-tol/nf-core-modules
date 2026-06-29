process LONGC_PAIRTOOLSPARSE2 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.1.3--py39h7a39fba_0' :
        'quay.io/biocontainers/pairtools:1.1.3--py39h7a39fba_0' }"

    input:
    tuple val(meta), path(bam), path(chromsizes)

    output:
    tuple val(meta), path("*.pairs.gz"), emit: pairs
    tuple val("${task.process}"), val('pairtools'), eval("pairtools --version | sed 's/.*pairtools.*version //'"), topic: versions, emit: versions_pairtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pairtools \\
        parse2 \\
        -c ${chromsizes} \\
        -o ${prefix}.pairs.gz \\
        ${args} \\
        ${bam}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.pairs.gz
    """
}
