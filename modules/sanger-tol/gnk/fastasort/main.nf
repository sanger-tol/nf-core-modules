process GNK_FASTASORT {
    tag "${meta.id}"

    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/308db3849b816a57ea2a8b1234c440948f57d428e5dffa6ec88cc46b23af1364/data'
        : 'community.wave.seqera.io/library/pip_gnk-fastasort:6f54e70aa0811c1a'}"

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
        --index ${index} \\
        --output ${prefix} \\
        ${args}
    """


    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
