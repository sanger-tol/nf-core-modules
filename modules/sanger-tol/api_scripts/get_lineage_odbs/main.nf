process API_SCRIPTS_GET_LINEAGE_ODBS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34a3fbdf8feee942542033e0c6961cd1efb0f4b2d4ec580982c653a74a73949d/data'
        : 'community.wave.seqera.io/library/requests:2.34.2--26550d4c3ba1cb75' }"

    input:
    tuple val(meta), path(fasta)
    path(odb_dir)
    val(odb_version)
    val(taxid)
    val(mode)
    val(specified_lineages)

    output:
    tuple val(meta), path("*.busco_odb.csv"), emit: csv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('get_odbs.py'), eval('get_odbs.py --version | cut -d" " -f2'), emit: versions_get_odb, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def odb_dir_path    = odb_dir ? "--odb_dir ${odb_dir}" : ""
    valid_mode          = mode ? "--mode ${mode}" : ""
    formatted_lineages  = specified_lineages && mode?.contains("specified") ? "--specified_lineages " + specified_lineages.tokenize(',').join(' ') : ""
    """
    get_odbs.py \\
        --taxid ${taxid} \\
        --odb_version ${odb_version} \\
        ${odb_dir_path} \\
        ${valid_mode} \\
        --file_out ${prefix}.busco_odb.csv \\
        ${args} \\
        ${formatted_lineages}
    """

    stub:
    def prefix          = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.busco_odb.csv
    """
}
