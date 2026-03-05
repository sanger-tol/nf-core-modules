process GNK_FASTASORT {
    tag "${meta.id}"

    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7570c98be77c9ca032c4b9abd1adfbc8a908db3d2662b9e0e8199f389250fff/data'
        : 'community.wave.seqera.io/library/pip_gnk_fastasort:cfccd490e04b325e'}"

    input:
    tuple val(meta), path(index)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('gnk_fastasort'), eval("fastasort -v | sed 's/fastasort. //g'"), topic: versions, emit: versions_gnkfastasort

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ""
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    fastasort \\
        --index ${fai} \\
        ${args}
    """


    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
