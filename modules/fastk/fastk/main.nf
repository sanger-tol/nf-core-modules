process FASTK_FASTK {
    tag "$meta.id"
    label 'process_low'

    def version = '0.001-c3'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/fastk:${version}"

    input:
    tuple val(meta), path(trimmedReads)
    val kmer

    output:
    tuple val(meta), path("*.hist"), emit: hist
    tuple val(meta), path("*.ktab*", hidden: true), emit: ktab, optional: true
    tuple val(meta), path("*.prof*", hidden: true), emit: prof, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ?: ''
    def cups = task.cups ?: ''

    """
    FastK -k$kmer -T1 -t$cups -N${prefix}_fk $trimmedReads
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $version
    END_VERSIONS
    """
}
