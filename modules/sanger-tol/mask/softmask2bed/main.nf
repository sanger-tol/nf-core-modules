process MASK_SOFTMASK2BED {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1'
        : 'biocontainers/seqtk:1.4--he4a0461_1'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('softmask2bed'), eval('echo 1.0'), emit: versions_softmask2bed, topic: versions
    tuple val("${task.process}"), val('seqtk'), eval("seqtk 2>&1 | sed -n 's/^Version: //p'"), emit: versions_seqtk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = fasta.extension == "gz" ? "<(zcat ${fasta})" : "${fasta}"
    // Use "seqtk seq" to convert all soft-masked regions (lowercase) to
    // hard-masked (N) regions, then use "seqtk gap" to report the gaps in BED
    // format. For this to work, I need to first replace all N in the input
    // sequences with a non-N (here A), so that only the masked regions are
    // reported in the output BED file.
    """
    sed '/^[^>]/s/N/A/g' ${input} \\
        | seqtk seq -x -n N \\
        | seqtk gap -l 1 - \\
        > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
