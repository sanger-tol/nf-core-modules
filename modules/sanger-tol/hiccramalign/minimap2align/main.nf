process HICCRAMALIGN_MINIMAP2ALIGN {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/711a743a18039d9ecfc6032fb1e124eea4a87bc2bf2ed47a1c5b8534f8193b1f/data' :
        'community.wave.seqera.io/library/minimap2_samtools_perl:65833242e5ffdf5c' }"

    input:
    tuple val(meta), path(cram), path(crai), val(chunkn), val(range), path(reference)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARNING: This module includes the filter_five_end.pl script as a module binary in
    // ${moduleDir}/resources/usr/bin/filter_five_end.pl. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or enable the option
    // nextflow.enable.moduleBinaries = true
    // in your nextflow.config file.
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def args5 = task.ext.args5 ?: ''
    def args6 = task.ext.args6 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools cat ${args1} -r "#:${range[0]}-${range[1]}" ${cram} |\\
        samtools fastq ${args2} - |\\
        minimap2 -t${task.cpus} ${args3} ${reference} - |\\
        filter_five_end.pl |\\
        samtools fixmate ${args4} - - |\\
        samtools view -h ${args5} |\\
        samtools sort ${args6} -@${task.cpus} -T ${prefix}_tmp -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
    END_VERSIONS
    """
}
