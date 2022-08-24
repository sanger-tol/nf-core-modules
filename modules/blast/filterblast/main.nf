
process BLAST_FILTERBLAST {
    tag "$meta.id"
    label 'process_medium'
    def version = '0.001-c2'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the filterblast process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/genealignment:${version}"


    input:
    tuple val( meta ), file( concat_blast_out )

    output:
    tuple val( meta ), file( "${meta.id}-${meta.type}-*.tsv")   , emit: final_tsv
    path "versions.yml"                                         , emit: versions

    script:
    def filt_percent = task.ext.args ?: 90.00
    """
    filter_blast.py $meta.id $meta.type $concat_blast_out $filt_percent
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_blast: $version
    END_VERSIONS
    """
}
