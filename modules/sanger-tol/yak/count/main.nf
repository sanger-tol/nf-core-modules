process YAK_COUNT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f4/f4f7c8cf06b8a0857b65315a6374dc76d86d18ce1aab04c4a8e8e3fb43ac9d07/data' :
        'community.wave.seqera.io/library/samtools_yak:c23b2980c0d714c3' }"

    input:
    tuple val(meta), path(reads, arity: '1..*')

    output:
    tuple val(meta), path("*.yak"), emit: yak
    tuple val("${task.process}"), val('yak'), eval("yak version"), topic: versions, emit: versions_yak

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cram_input = (reads[0].getExtension() in ["bam", "cram"])
    def fasta_command = cram_input ? "samtools cat ${reads} | samtools fasta" : "zcat ${reads}"
    """
    yak \\
        count \\
        $args \\
        -t${task.cpus} \\
        -o ${prefix}.yak \\
        <(${fasta_command}) <(${fasta_command})
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.yak
    """
}
