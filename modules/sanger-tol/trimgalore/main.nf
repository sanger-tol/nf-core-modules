process TRIMGALORE {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_low_memory'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92076024c2fd5addd785ce04fef5f754fb744c668efb457ff19f20150a405999/data'
        : 'community.wave.seqera.io/library/trim-galore_samtools_gzip:270d120a8d462d4a'}"

    input:
    tuple val(meta), path(cram)

    output:
    tuple val(meta), path("*{trimmed,val}.cram")                       , emit: cram 
    tuple val(meta), path("*report.txt")                               , emit: log     , optional: true
    tuple val(meta), path("*.html")                                    , emit: html    , optional: true
    tuple val(meta), path("*.zip")                                     , emit: zip     , optional: true
    tuple val("${task.process}"), val("trimgalore"), eval('trim_galore --version | grep -Eo "[0-9]+(\\.[0-9]+)+"'), topic: versions, emit: versions_trimgalore

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/CHANGELOG.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (cores < 1) {
            cores = 1
        }
        if (cores > 8) {
            cores = 8
        }
    }

    // Convert from cram to unaligned bam and pass to trim_galore 
    // Will likely not work for single-end RNAseq datasets - not in scope of current workflow

    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -h -@ ${task.cpus} -O BAM ${cram} -o ${prefix}.bam 

    trim_galore \\
        ${args} \\
        --cores ${cores} \\
        --paired \\
        --output-format ubam \\
        ${prefix}.bam

    samtools view -h -@ ${task.cpus} -C ${prefix}_val.bam \\
        -o ${prefix}_val.cram

    rm ${prefix}*.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_command = "echo '' | gzip > ${prefix}_val.cram ;"
    output_command += "touch ${prefix}_1.fastq.gz_trimming_report.txt ;"
    output_command += "touch ${prefix}_2.fastq.gz_trimming_report.txt"
    """
    ${output_command}
    """
}
