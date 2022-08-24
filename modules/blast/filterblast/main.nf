
process BLAST_FILTERBLAST {
    tag "$meta.id"
    label 'process_medium'
    def version = '0.001-c2'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the filterblast process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/genealignment:${version}"


    input:
    tuple val(meta), path(blastout)

    output:
    tuple val( meta ), path( "${meta.id}-${meta.type}-*.tsv")   , emit: final_tsv
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def filt_percent = task.ext.args ?: 90.00
    """
    filter_blast.py $meta.id $meta.type $blastout $filt_percent
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_blast: $version
    END_VERSIONS
    """
}
