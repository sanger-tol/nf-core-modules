process NCBIDATASETS_GET_ODB {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::requests=2.26.0"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/requests:2.26.0'
        : 'community.wave.seqera.io/library/requests:2.34.2--26550d4c3ba1cb75' }"

    input:
    tuple val(meta), path(ncbi_summary)
    path lineage_tax_ids
    val all_ancestral_lineages

    output:
    tuple val(meta), path("*.busco_odb.csv"), emit: csv
    tuple val("${task.process}"), val('python'), eval('Python --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('get_odb.py'), eval('get_odb.py --version | cut -d" " -f2'), emit: versions_get_odb, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    get_odb.py \\
        --ncbi_summary_json ${ncbi_summary} \\
        --lineage_tax_ids ${lineage_tax_ids} \\
        --all_ancestral_lineages ${all_ancestral_lineages} \\
        ${args} \\
        ${prefix}.busco_odb.csv
    """

    stub:
    """
    touch ${prefix}.busco_odb.csv
    """
}
