process AUTOFILTER_SANITISEFASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(input_fasta)

    output:
    tuple val(meta), path("*_filtered.fasta"),              emit: fasta
    tuple val(meta), path("fasta_sanitation.json"),         emit: sanitation_log
    tuple val(meta), path("fasta_length_filtering.json"),   emit: length_filtering_log

    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('sanitise_fasta_headers'), eval('sanitise_fasta_headers.py --version | cut -d" " -f2'), emit: versions_sanitise,  topic: versions
    tuple val("${task.process}"), val('filter_fasta_by_length'), eval('filter_fasta_by_length.py --version | cut -d" " -f2'), emit: versions_length,    topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    def args        = task.ext.args     ?: ''
    def args2       = task.ext.args2    ?: ''
    """
    sanitise_fasta_headers.py \\
        ${input_fasta} \\
        --log_file fasta_sanitation.json \\
        --max_detailed_changes 0 > ${prefix}_shortened.fasta \\
        ${args}

    filter_fasta_by_length.py \\
        ${args2} \\
        ${prefix}_shortened.fasta \\
        --log_file fasta_length_filtering.json > ${prefix}_filtered.fasta
    """

    stub:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_filtered.fasta
    touch fasta_sanitation.json
    touch fasta_length_filtering.json
    """
}
