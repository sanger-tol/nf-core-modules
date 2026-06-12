process PRETEXT_PRETEXT2ASM {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/sanger-tol/agp-tpf-utils:1.3.4-c2'

    input:
    tuple val(meta), path(mappedfasta)
    path(pretextagp)


    output:
    tuple val(meta), path("${meta.id}_corrected*.agp"), emit: correctedagp
    // tola-agp-tpf-utils has no CLI --version; pin to container tag (bump when `container` changes)
    tuple val("${task.process}"), val('tola-agp-tpf-utils'), val('1.3.4'), topic: versions, emit: versions_agp_tpf_utils

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PRETEXT_PRETEXT2ASM module does not support Conda. Please use Docker / Singularity instead."
    }
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    pretext-to-asm -a ${mappedfasta} -p ${pretextagp} -o ${prefix}_corrected.agp ${args}
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PRETEXT_PRETEXT2ASM module does not support Conda. Please use Docker / Singularity instead."
    }
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_corrected.agp
    """
}
