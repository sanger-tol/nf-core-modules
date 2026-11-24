process CRAMALIGN_STARALIGN {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/faf2d70727cd58d49615fd0679019bb7029e6fc04ec8eb2c42d510d24b0af19a/data'
        : 'community.wave.seqera.io/library/htslib_samtools_star_gawk:2a8b56d2f32c7260'}"

    input:
    tuple val(meta),  path(cram),  path(crai), val(rglines)
    tuple val(meta2), path(index), path(assembly)
    tuple val(chunkn), val(range)

    output:
    tuple val(meta), path("*.star.bam"), emit: bam
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed -e "s/STAR_//g"'), emit: versions_star, topic: versions
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -1 | sed -e "s/samtools //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${cram}.${chunkn}.${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ""
    def args4 = task.ext.args4 ?: ''
    // Prepare read group arguments if rglines are found, else, empty string
    def rg_arg = rglines ? '-C ' + rglines.collect { line ->
            // Add SM when not present to avoid errors from downstream tool (e.g. variant callers)
            def l = line.contains("SM:") ? line : "${line}\tSM:${meta.id}"
            "-H '${l.replaceAll("\t", "\\\\t")}'"
        }.join(' ')
        : ''
    """

    samtools view -b -h -o unaligned.bam ${cram}
    STAR \\
        --genomeDir $index \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $args \\
        --readFilesCommand "samtools view -h" \\
        --readFilesIn unaligned.bam

    samtools view -h ${prefix}.Aligned.out.bam |
        insert_cram_pg_header.awk -v pgfile="${prefix}_cram_pg.tmp" |\\
        samtools fixmate ${args2} - - |\\
        samtools view -h ${args3} |\\
        samtools sort ${args4} -@${task.cpus} -T ${prefix}_tmp -o ${prefix}.star.bam -

    rm unaligned.bam
    rm ${prefix}.Aligned.out.bam
    """

    stub:
    def prefix  = task.ext.prefix ?: "${cram}.${chunkn}.${meta.id}"
    """
    touch ${prefix}.star.bam
    """
}
