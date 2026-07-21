process AUTOFILTER_AUTOFILTER {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.14' :
        'quay.io/biocontainers/python:3.14' }"

    input:
    tuple val(meta),        path(fai)
    tuple val(tiara_meta),  path(tiara_txt)
    tuple val(fcs_meta),    path(fcs_csv)
    val taxid
    path ncbi_rankedlineage

    output:
    tuple val(meta), path("*_autofiltered.txt"),       emit: keep_scaffs
    tuple val(meta), path("*_removed_scaffolds.txt"),  emit: remove_scaffolds
    tuple val(meta), path("*_ABNORMAL_CHECK.csv"),     emit: fcs_tiara_summary
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('autofilter'), eval('autofilter.py --version | cut -d" " -f2'), emit: versions_autofilter, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ""
    """
    autofilter.py \\
        $fai \\
        --taxid $taxid \\
        --tiara $tiara_txt \\
        --fcsgx_sum $fcs_csv \\
        --prefix $prefix \\
        --ncbi_rankedlineage_path $ncbi_rankedlineage \\
        ${args}
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    touch ${prefix}_autofiltered.txt
    touch ${prefix}_removed_scaffolds.txt
    touch ${prefix}_ABNORMAL_CHECK.csv
    """
}
