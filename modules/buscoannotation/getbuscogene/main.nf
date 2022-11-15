process BUSCOANNOTATION_GETBUSCOGENE {
    tag "$meta.id"
    label 'process_medium'

    def version = '0.001-c5'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the buscogeneannotation process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/buscogeneannotation:${version}"
    
    input:
    tuple val(meta), path(full_table_tsv)

    output:
    tuple val(meta), path("*_busco.csv")  , emit: busco_csv
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_busco_gene.py \\
        $args \\
        $full_table_tsv \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        buscoannotation: $version
    END_VERSIONS
    """
}