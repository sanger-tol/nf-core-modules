process AUTOFILTER_CONTAMINATIONCHECK {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta),    path(reference)
    tuple val(meta2),   path(check_csv)

    output:
    tuple val(meta), path("alarm_indicator_file.txt"),              emit: alarm_file
    tuple val(meta), path("autofiltering_done_indicator_file.txt"), emit: indicator_file
    tuple val(meta), path("*raw_report.txt"),                       emit: raw_report
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('abnormal_contamination_check'), eval('abnormal_contamination_check.py --version | cut -d" " -f2'), emit: versions_autofilter, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ""
    """
    abnormal_contamination_check.py \\
        ${reference} \\
        ${check_csv} \\
        --out_prefix ${prefix}

    # The indicator file is used in Sanger-Tol to allow for other processes
    # to begin once generated. This allows us to speed up the overall flow of the
    # Tol-engine
    echo "AUTOFILTER COMPLETE" > autofiltering_done_indicator_file.txt
    """

    stub:
    """
    touch alarm_indicator_file.txt
    touch autofiltering_done_indicator_file.txt
    touch raw_report.txt
    """
}
