process AUTOFILTER_FILTERFASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.14' :
        'quay.io/biocontainers/python:3.14' }"

    input:
    tuple val(meta), path(fai)

    output:
    tuple val(meta), path("*_mapping.txt"),         emit: mapping_file
    tuple val(meta), path("*_new_scaffolds.txt"),   emit: new_scaffolds
    tuple val(meta), path("*_stats.json"),          emit: stats

    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('filter_fasta'), eval('filter_fasta.py --version | cut -d" " -f2'), emit: versions_filter_fasta,  topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    def args        = task.ext.args     ?: ''
    def args2       = task.ext.args2    ?: ''
    """
    filter_fasta.py \\
        ${fai} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_mapping.txt
    touch ${prefix}_new_scaffolds.txt
    touch ${prefix}_stats.json
    """
}
