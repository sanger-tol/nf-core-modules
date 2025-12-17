process CRAMALIGN_MINIMAP2ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65858e733832166824cfd05291fc456bdf219b02baa3944c2c92efad86a6ee7f/data' :
        'community.wave.seqera.io/library/htslib_minimap2_samtools_gawk_perl:6729620c63652154' }"

    input:
    tuple val(meta), val(rglines), path(cram), path(crai), val(chunkn), val(range), path(reference)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: '-t' // copy RG, BC and QT tags to the FASTQ header line
    def args3 = task.ext.args3 ?: ''
    def args5 = task.ext.args5 ?: ''
    def prefix = task.ext.prefix ?: "${cram}.${chunkn}.${meta.id}"
    // Allow full region if range not provided
    def (samtools_cat_cmd, cram_input) = range
        ? [ "samtools cat ${args1} -r '#:${range[0]}-${range[1]}' ${cram} |", "-" ]
        : [ '', cram ]
    def post_filter = task.ext.args4 ? "samtools view -h ${task.ext.args4} |" : ''
    def rg_arg = rglines ? '-y ' + rglines.collect { line ->
            // Add SM when not present to avoid errors from downstream tool (e.g. variant callers)
            def l = line.contains("SM:") ? line
                : meta.sample ? "${line}\tSM:${meta.sample}"
                : "${line}\tSM:${meta.id}"
            "-R '${l.replaceAll("\t", "\\\\t")}'"
        }.join(' ')
        : ''
    """
    ${samtools_cat_cmd} \\
        samtools fastq ${args2} ${cram_input} |  \\
        minimap2 -t${task.cpus} ${args3} ${reference} ${rg_arg} - | \\
        ${post_filter} \\
        samtools sort ${args5} -@${task.cpus} -T ${prefix}_sort_tmp -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${cram}.${chunkn}.${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
    END_VERSIONS
    """
}
