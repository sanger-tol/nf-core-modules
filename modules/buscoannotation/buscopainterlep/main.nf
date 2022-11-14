process BUSCOANNOTATION_BUSCOPAINTERLEP {
    tag "$meta.id"
    label 'process_medium'

    def version = '0.001-c4'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the busco painter process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/buscogeneannotation:${version}"

    input:
    tuple val(meta), path(full_table_tsv)
    path(merian_tsv)

    output:
    tuple val( meta ), path( "*_complete_location.tsv" )    , emit: complete_location_tsv
    tuple val( meta ), path( "*_duplicated_location.tsv" )  , emit: duplicated_location_tsv
    tuple val( meta ), path( "*_summary.tsv" )              , emit: summary_tsv
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lep_buscopainter.py \\
        -r $merian_tsv \\
        -q $full_table_tsv \\
        -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        buscogeneannotation: ${version}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_complete_location.tsv
    touch ${prefix}_duplicated_location.tsv
    touch ${prefix}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        buscogeneannotation: ${version}
    END_VERSIONS
    """
}
