process PRETEXT_MAKEPAIRS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(alignment), path(outlog)

    output:
    tuple val(meta), path("*_alignment.pairs.gz"), emit: pairs
    // zippypretext make_header / makepairs scripts have no CLI --version; pin to 1.0.0
    tuple val("${task.process}"), val('pretext_makepairs'), val('1.0.0'), topic: versions, emit: versions_pretext_makepairs

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    (
        echo '## pairs format v1.0'
        grep PRE_C_SIZE ${outlog} \\
            | tr -s ' ' '\\t' \\
            | cut -f2,3 \\
            | while IFS=\$'\\t' read -r chr size; do
                printf '#chromsize:\\t%s\\t%s\\n' "\$chr" "\$size"
            done
        printf '#columns:\\treadID\\tchr1\\tpos1\\tchr2\\tpos2\\tstrand1\\tstrand2\\n'
        cut -f2,3,6,7 ${alignment} \\
            | while IFS=\$'\\t' read -r c1 p1 c2 p2; do
                printf '.\\t%s\\t%s\\t%s\\t%s\\t.\\t.\\n' "\$c1" "\$p1" "\$c2" "\$p2"
            done
    ) | gzip > ${prefix}_alignment.pairs.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_alignment.pairs.gz
    """
}
