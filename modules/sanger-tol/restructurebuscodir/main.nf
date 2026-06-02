process RESTRUCTUREBUSCODIR {
    tag "${meta.id}_${lineage}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), val(lineage), path(batch_summary), path(short_summary_txt), path(short_summary_json), path(full_table), path(missing_busco_list), path(seq_dir)

    output:
    tuple val(meta), path("${lineage}"), emit: clean_busco_dir
    tuple val("${task.process}"), val('restructurebuscodir'), eval('echo 1.0'), topic: versions, emit: versions_restructurebuscodir

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${lineage}

    [ -e "${batch_summary}" ]      && ln -s ../${batch_summary}      ${lineage}/${prefix}.short_summary.tsv
    [ -e "${short_summary_txt}" ]  && ln -s ../${short_summary_txt}  ${lineage}/${prefix}.short_summary.txt
    [ -e "${short_summary_json}" ] && ln -s ../${short_summary_json} ${lineage}/${prefix}.short_summary.json

    [ -e "${full_table}" ]         && gzip -c < ${full_table}         > ${lineage}/${prefix}.full_table.tsv.gz
    [ -e "${missing_busco_list}" ] && gzip -c < ${missing_busco_list} > ${lineage}/${prefix}.missing_busco_list.tsv.gz

    if [ -e "${seq_dir}/single_copy_busco_sequences" ]
    then
        tar czf ${lineage}/${prefix}.single_copy_busco_sequences.tar.gz -C ${seq_dir} single_copy_busco_sequences
        tar czf ${lineage}/${prefix}.multi_copy_busco_sequences.tar.gz  -C ${seq_dir} multi_copy_busco_sequences
        tar czf ${lineage}/${prefix}.fragmented_busco_sequences.tar.gz  -C ${seq_dir} fragmented_busco_sequences
    else
        # Busco was run in --tar mode, the sequences are already compressed
        ln -s ../${seq_dir}/single_copy_busco_sequences.tar.gz ${lineage}/${prefix}.single_copy_busco_sequences.tar.gz
        ln -s ../${seq_dir}/multi_copy_busco_sequences.tar.gz ${lineage}/${prefix}.multi_copy_busco_sequences.tar.gz
        ln -s ../${seq_dir}/fragmented_busco_sequences.tar.gz ${lineage}/${prefix}.fragmented_busco_sequences.tar.gz
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${lineage}
    touch ${lineage}/${prefix}.short_summary.tsv
    touch ${lineage}/${prefix}.short_summary.txt
    touch ${lineage}/${prefix}.short_summary.json
    touch ${lineage}/${prefix}.full_table.tsv
    touch ${lineage}/${prefix}.missing_busco_list.tsv
    gzip ${lineage}/${prefix}.full_table.tsv
    gzip ${lineage}/${prefix}.missing_busco_list.tsv
    tar -czf ${lineage}/${prefix}.single_copy_busco_sequences.tar.gz -T /dev/null
    tar -czf ${lineage}/${prefix}.multi_copy_busco_sequences.tar.gz  -T /dev/null
    tar -czf ${lineage}/${prefix}.fragmented_busco_sequences.tar.gz  -T /dev/null
    """
}
