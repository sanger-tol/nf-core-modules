process HICCRAMALIGN_BWAMEM2ALIGN {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/41a128095da89d592889899dceed75ba0252bf9471eebee2db6ae0eb5474738a/data' :
        'community.wave.seqera.io/library/bwa-mem2_samtools:34483fa65c9842a5' }"

    input:
    tuple val(meta), path(cram), path(crai), val(chunkn), val(range), path(index), path(reference)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def args5 = task.ext.args5 ?: ''
    def args6 = task.ext.args6 ?: ''
    def prefix  = task.ext.prefix ?: "${cram}.${chunkn}.${meta.id}"
    // Please be aware one of the tools here required mem = 28 * reference size!!!
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    samtools cat ${args1} -r "#:${range[0]}-${range[1]}" ${cram} |\\
        samtools fastq ${args2} - |\\
        bwa-mem2 mem ${args3} -t ${task.cpus} \${INDEX} - |\\
        samtools fixmate ${args4} - - |\\
        samtools view -h ${args5} |\\
        samtools sort ${args6} -@${task.cpus} -T ${prefix}_tmp -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${cram}.${chunkn}.${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
