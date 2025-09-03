process BLOBTK_DEPTH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/genomehubs/blobtk:0.6.5"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*.regions.bed.gz') , emit: bed
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    blobtk depth \\
        -b ${bam} \\
        $args \\
        -O ${prefix}.regions.bed.gz \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtk: \$(blobtk --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION     = "0.6.5"
    """
    echo "" | gzip > ${prefix}.regions.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtk: $VERSION
    END_VERSIONS
    """
}
