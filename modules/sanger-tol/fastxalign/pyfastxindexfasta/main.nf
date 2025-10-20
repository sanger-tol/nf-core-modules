process FASTXALIGN_PYFASTXINDEXFASTA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65858e733832166824cfd05291fc456bdf219b02baa3944c2c92efad86a6ee7f/data' :
        'community.wave.seqera.io/library/htslib_minimap2_samtools_gawk_perl:6729620c63652154' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fxi"), stdout, emit: index
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARNING: This module includes the slice_fasta.py script as a module binary in
    // ${moduleDir}/resources/usr/bin/slice_fasta.py. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or set the option
    // nextflow.enable.moduleBinaries = true
    // in your nextflow.config file.
    def args       = task.ext.args  ?: ''
    """
    slice_fasta.py index ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slice_fasta.py: \$(slice_fasta.py --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.fxi
    ## output dummy count to stdout
    echo -n 100

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slice_fasta.py: \$(slice_fasta.py --version)
    END_VERSIONS
    """
}
