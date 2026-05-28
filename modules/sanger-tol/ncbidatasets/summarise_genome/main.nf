process NCBIDATASETS_SUMMARISE_GENOME {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "docker.io/biocontainers/ncbi-datasets-cli:16.22.1_cv1"

    errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json"), emit: summary
    tuple val("${task.process}"), val('python'), eval('Python --version | sed "s/Python //"'), emit: versions_python, topic: versions
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
