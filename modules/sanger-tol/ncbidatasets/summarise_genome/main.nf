process NCBIDATASETS_SUMMARISE_GENOME {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65bfaf0aeeed8b31e36bd8effe70d1afb8fab6768e015fe25834dfa3449a7fb6/data'
        : 'community.wave.seqera.io/library/ncbi-datasets-cli_python:9f4299728bfbaaa1' }"


    errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json"), emit: summary
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val("datasets"), eval("datasets --version | sed 's/^.*datasets version. //'"), emit: versions_datasets, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    datasets \\
        summary \\
        genome \\
        accession \\
        ${meta.id} \\
        ${args} \\
        > ${prefix}.json

    validate_datasets_json.py ${prefix}.json
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json
    """
}
