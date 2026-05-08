process FINDTELOMERE {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/sanger-tol/telomere:0.0.2-c1'

    input:
    tuple val(meta), path(reference), val(telomereseq)
    val split_windows

    output:
    tuple val(meta), path("*.telomere"), emit: telomere
    tuple val(meta), path("*.fwd.telomere.bed"), emit: telomere_bed_fwd, optional: true
    tuple val(meta), path("*.rev.telomere.bed"), emit: telomere_bed_rev, optional: true
    tuple val(meta), path("*.windows"), emit: windows, optional: true
    tuple val(meta), path("*.fwd.windows"), emit: windows_fwd, optional: true
    tuple val(meta), path("*.rev.windows"), emit: windows_rev, optional: true
    tuple val("${task.process}"), val('java'), eval("java -version 2>&1 | head -n 1 | cut -d '\"' -f2"), topic: versions, emit: versions_java
    tuple val("${task.process}"), val('find_telomere'), eval("find_telomere 2>&1 | head -n 1 | sed 's/^[[:space:]]*//'"), topic: versions, emit: versions_find_telomere

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FINDTELOMERE module does not support Conda. Please use Docker / Singularity instead."
    }
    def args = task.ext.args ?: '{print \$1"\\t"\$(NF-4)"\\t"\$(NF-3)"\\t"\$(NF-2)"\\t"\$(NF-1)"\\t"\$NF}'
    def args2 = task.ext.args2 ?: '-Xmx4096M -Xss999M'
    def args3 = task.ext.args3 ?: '99.9 0.1'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def split_opt = split_windows ? '--split ' : ''
    def stdout_redirect = split_windows ? '' : "> ${prefix}.windows"
    """
    find_telomere $reference $telomereseq | awk '${args}' - > ${prefix}.telomere

    java \\
        ${args2} \\
        -cp /opt/telomere/telomere.jar \\
        FindTelomereWindows \\
        ${split_opt}${prefix}.telomere \\
        ${args3} \\
        ${stdout_redirect}
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FINDTELOMERE module does not support Conda. Please use Docker / Singularity instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "stub\\n" > ${prefix}.telomere
    printf "stub\\n" > ${prefix}.fwd.telomere.bed
    printf "stub\\n" > ${prefix}.rev.telomere.bed
    if ${split_windows}; then
        printf "stub\\n" > ${prefix}.fwd.windows
        printf "stub\\n" > ${prefix}.rev.windows
    else
        printf "stub\\n" > ${prefix}.windows
    fi
    """

}
