process PRETEXT_JUICERC {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/sanger-tol/juicerc:1.2-c1'

    input:
    tuple val(meta), path(hicmap)
    path agpfile
    path idxfile

    output:
    tuple val(meta), path("*_alignment_sorted.txt"), emit: alignment
    path("juicercout.log"), emit: outlog
    tuple val("${task.process}"), val('juicer_pre'), eval('juicer pre --version'), topic: versions, emit: versions_juicer_pre
    // juicerc container has no CLI --version; pin to container tag (bump when `container` changes)
    tuple val("${task.process}"), val('juicerc'), val('1.2'), topic: versions, emit: versions_juicerc

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PRETEXT_JUICERC module does not support Conda. Please use Docker / Singularity instead."
    }
    def prefix   = task.ext.prefix   ?: "${meta.id}"
    def sort_mem = (task.memory.toGiga() * 0.9).intValue()

    """
    juicer pre -q 0 ${hicmap} ${agpfile} ${idxfile} 2>juicercout.log \\
        | LC_ALL=C sort -k2,2d -k6,6d -T . --parallel=${task.cpus} -S${sort_mem}G \\
        | awk '\$3>=0 && \$7>=0' \\
        > ${prefix}_alignment_sorted.txt
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PRETEXT_JUICERC module does not support Conda. Please use Docker / Singularity instead."
    }
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_alignment_sorted.txt
    touch juicercout.log
    """
}
