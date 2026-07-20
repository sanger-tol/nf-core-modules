process AUTOFILTER_AUTOFILTER {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.14' :
        'quay.io/biocontainers/python:3.14' }"

    input:
    tuple val(meta),        path(reference)
    tuple val(tiara_meta),  path(tiara_txt)
    tuple val(fcs_meta),    path(fcs_csv)
    val taxid
    path ncbi_rankedlineage

    output:
    tuple val(meta), path("*autofiltered.fasta"),                       emit: decontaminated_assembly, optional: true
    tuple val(meta), path("*ABNORMAL_CHECK.csv"),                       emit: fcs_tiara_summary
    tuple val(meta), path("assembly_filtering_removed_sequences.txt"),  emit: removed_seqs
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('autofilter'), eval('autofilter.py --version | cut -d" " -f2'), emit: versions_autofilter, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ""
    """
    autofilter.py \\
        $reference \\
        --taxid $taxid \\
        --tiara $tiara_txt \\
        --fcsgx_sum $fcs_csv \\
        --out_prefix $prefix \\
        --ncbi_rankedlineage_path $ncbi_rankedlineage \\
        ${args}
    """

    stub:
    """
    touch autofiltered.fasta
    touch ABNORMAL_CHECK.csv
    touch assembly_filtering_removed_sequences.txt
    """
}
