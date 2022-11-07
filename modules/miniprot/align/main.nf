
process MINIPROT_ALIGN {

    tag "$meta.id"
    label 'process_medium'

    def version = '0.5-c2'
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the miniprot process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/miniprot:${version}"

    input:
    tuple val(meta), path(ref)
    path pep                     /* Can be fasta format for index */
    val paf_format
    val gff_format
    val gtf_format

    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.gff"), optional: true, emit: gff
    tuple val(meta), path("*.gtf"), optional: true, emit: gtf
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paf_out = paf_format ? "> ${prefix}.paf" : ''
    def gff_out = gff_format ? "> ${prefix}.gff" : ''
    def gtf_out = gtf_format ? "> ${prefix}.gtf" : ''
    """
    miniprot \\
        $args \\
        -t $task.cpus \\
        $ref \\
        $pep \\
        $paf_out \\
        $gff_out \\
        $gtf_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        miniprot: $version
    END_VERSIONS
    """
}

