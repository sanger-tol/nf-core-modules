process PRETEXT_MAKEPAIRS {
    tag "$meta.id"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'quay.io/biocontainers/ubuntu:20.04'}"

    input:
    tuple val(meta), path(alignment), path(outlog)

    output:
    tuple val(meta), path("*_alignment.pairs.gz"), emit: pairs
    // zippypretext make_header / makepairs scripts have no CLI --version; pin to 1.0.0
    tuple val("${task.process}"), val('pretext_makepairs'), val('1.0.0'), topic: versions, emit: versions_pretext_makepairs

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PRETEXT_MAKEPAIRS module does not support Conda. Please use Docker / Singularity instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    (
        grep PRE_C_SIZE ${outlog} \\
            | awk '{print \$2"\\t"\$3}' \\
            | awk 'BEGIN{print "## pairs format v1.0"} {print "#chromsize:\\t"\$1"\\t"\$2} END{print "#columns:\\treadID\\tchr1\\tpos1\\tchr2\\tpos2\\tstrand1\\tstrand2"}'
        awk '{print ".\\t"\$2"\\t"\$3"\\t"\$6"\\t"\$7"\\t.\\t."}' ${alignment}
    ) | gzip > ${prefix}_alignment.pairs.gz
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PRETEXT_MAKEPAIRS module does not support Conda. Please use Docker / Singularity instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_alignment.pairs.gz
    """
}
