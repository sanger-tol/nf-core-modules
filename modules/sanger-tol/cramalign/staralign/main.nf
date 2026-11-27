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
    tuple val(meta), path("*.bam"), emit: bam
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
    def args5 = task.ext.args4 ?: ''
    """
    samtools view --no-PG -H ${cram} | grep ^@PG > ${prefix}_cram_pg.tmp

    samtools cat ${args} -r "#:${range[0]}-${range[1]}" ${cram} |\\
    samtools view -b -h -o ${cram.getBaseName()}.bam -

    STAR \\
        --genomeDir $index \\
        --runThreadN $task.cpus \\
        --outSAMtype BAM Unsorted \\
        --readFilesType SAM PE \\
        --outStd BAM_Unsorted \\
        ${args2} \\
        --readFilesCommand "samtools view -h" \\
        --readFilesIn ${cram.getBaseName()}.bam |\\
    insert_cram_pg_header.awk -v pgfile="${prefix}_cram_pg.tmp" |\\
    samtools fixmate ${args3} - - |\\
    samtools view -h ${args4} |\\
    samtools sort ${args5} -@${task.cpus} -T ${prefix}_tmp -o ${prefix}.bam -

    rm ${cram.getBaseName()}.bam
    """

    stub:
    def prefix  = task.ext.prefix ?: "${cram}.${chunkn}.${meta.id}"
    """
    touch ${prefix}.bam
    """
}
