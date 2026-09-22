process BGZIPTABIX {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/86/863ca0dbbba30c8367fa4fbd3fa3a84393532fb7b300a5c5c2e70f0dfc475bbf/data'
        : 'community.wave.seqera.io/library/htslib_xz:32f2772a564b3cd2'}"

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
    extension ?= input.extension
    """
    filter_compress () {
        if [[ -z "${column_numbers}" ]]
        then
            # No column / line filtering
            bgzip --threads ${task.cpus} --index ${args} --output ${prefix}.${extension}.gz
        else
            # Column / line filtering
            cut -f${column_numbers} | tail -n+${header_lines + 1} | bgzip --threads ${task.cpus} --index ${args} --output ${prefix}.${extension}.gz
        fi
    }

    FILE_TYPE=\$(htsfile ${input})

    DECOMPRESS=()
    NEED_COMPRESS=1

    case "\$FILE_TYPE" in
        *BGZF-compressed*)
            if [[ -z "${column_numbers}" ]]
            then
                ln -s "${input}" "${prefix}.${extension}.gz"
                bgzip --threads ${task.cpus} --reindex ${args} "${prefix}.${extension}.gz"
                NEED_COMPRESS=0
            else
                DECOMPRESS=(bgzip -d -c -@ "${task.cpus}")
            fi
            ;;
        *gzip-compressed*)
            DECOMPRESS=(bgzip -d -c -@ "${task.cpus}")
            ;;
        *bzip2-compressed*)
            DECOMPRESS=(bzcat)
            ;;
        *XZ-compressed*)
            DECOMPRESS=(xzcat)
            ;;
        *)
            ;;
    esac

    if ((NEED_COMPRESS))
    then
        if ((\${#DECOMPRESS[@]}))
        then
            filter_compress < <("\${DECOMPRESS[@]}" ${input})
        else
            filter_compress < ${input}
        fi
    fi

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
