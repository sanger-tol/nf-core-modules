process BGZIPTABIX {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data'
        : 'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa'}"

    input:
    tuple val(meta), path(input)
    tuple val(column_numbers), val(header_lines), val(sort_columns), val(extension)
    tuple val(tabix_tbi), val(tabix_csi), val(max_seq_length)

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
    extension ?= input.extension
    def output = "${prefix}.${extension}.gz"
    def filter_cut = column_numbers ? "cut -f${column_numbers} | " : ""
    def filter_tail = header_lines ? "tail -n+${header_lines + 1} | " : ""
    def filter_sort = sort_columns ? "sort ${sort_columns} | " : ""
    def filter = filter_cut + filter_tail + filter_sort
    def compress = "bgzip --threads ${task.cpus} --index ${args} --output ${output}"
    def do_tbi = tabix_tbi ? "1" : ""
    def do_csi = tabix_csi ? "1" : ""
    // default to not checking the size
    def msl = max_seq_length ?: 0
    """
    # Copied from the nf-core samtools/bgzip module

    FILE_TYPE=\$(htsfile ${input})
    case "\$FILE_TYPE" in
        *BGZF-compressed*)
            if [[ "${filter}" != "" ]]
            then
              zcat ${input} | ${filter} ${compress}
            else
              # Do nothing or just rename if the file was already compressed
              [[ "\$(basename ${input})" != "\$(basename ${output})" ]] && ln -s ${input} ${output}
            fi;;
        *gzip-compressed*)
            [[ "\$(basename ${input})" == "\$(basename ${output})" ]] && echo "Filename collision (\$basename ${input})" && exit 1
            zcat  ${input} | ${filter} ${compress};;
        *bzip2-compressed*)
            bzcat ${input} | ${filter} ${compress};;
        *XZ-compressed*)
            xzcat ${input} | ${filter} ${compress};;
        *)
            < ${input} ${filter} ${compress};;
    esac

    if [[ "${do_tbi}" != "" ]]
    then
        [[ ${msl} -lt \$(( 2 ** 29 )) ]] && tabix --threads ${task.cpus} ${args2} ${prefix}.${extension}.gz
    fi
    if [[ "${do_csi}" != "" ]]
    then
        [[ ${msl} -lt \$(( 2 ** 32 )) ]] && tabix --threads ${task.cpus} --csi ${args2} ${prefix}.${extension}.gz
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension ?= input.extension
    """
    echo "" | gzip > ${prefix}.${extension}.gz
    touch ${prefix}.${extension}.gz.gzi
    touch ${prefix}.${extension}.gz.tbi
    touch ${prefix}.${extension}.gz.csi
    """
}
