process BGZIPTABIX {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(input), val(max_seq_length)
    tuple val(column_numbers), val(header_lines), val(extension)

    output:
    tuple val(meta), path("*.gz"), path("*.gzi"), emit: gz_index
    tuple val(meta), path("*.tbi"), emit: tbi, optional: true
    tuple val(meta), path("*.csi"), emit: csi, optional: true
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_bgzip
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_data = column_numbers ? "<(cut -f${column_numbers} ${input} | tail -n+${header_lines + 1})" : input
    extension ?= input.extension
    """
    bgzip --threads ${task.cpus} --index ${args} ${input_data} --output ${prefix}.${extension}.gz
    [[ ${max_seq_length} -lt \$(( 2 ** 29 )) ]] && tabix --threads ${task.cpus} ${args2} ${prefix}.${extension}.gz
    [[ ${max_seq_length} -lt \$(( 2 ** 32 )) ]] && tabix --threads ${task.cpus} --csi ${args2} ${prefix}.${extension}.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension ?= input.extension
    """
    echo "" | bgzip > ${prefix}.${extension}.gz
    touch ${prefix}.${extension}.gz.gzi
    touch ${prefix}.${extension}.gz.tbi
    touch ${prefix}.${extension}.gz.csi
    """
}
