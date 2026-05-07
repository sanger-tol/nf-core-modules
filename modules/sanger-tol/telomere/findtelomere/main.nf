process FINDTELOMERE {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/sanger-tol/telomere:0.0.2-c1'

    input:
    tuple val(meta), path(reference), val(telomereseq)
    val split_windows

    output:
    tuple val( meta ), path( "*.telomere" ) , emit: telomere
    tuple val( meta ), path( "*.fwd.telomere.bed" )  , optional: true, emit: telomere_bed_fwd
    tuple val( meta ), path( "*.rev.telomere.bed" )  , optional: true, emit: telomere_bed_rev
    tuple val( meta ), path( "*.windows" )     , optional: true, emit: windows
    tuple val( meta ), path( "*.fwd.windows" ) , optional: true, emit: windows_fwd
    tuple val( meta ), path( "*.rev.windows" ) , optional: true, emit: windows_rev
    tuple val("${task.process}"), val('java'), eval("java -version 2>&1 | head -n 1 | cut -d '\"' -f2"), topic: versions, emit: versions_java
    /*
     * Avoid calling `find_telomere --version`: some builds treat `--version` as an input filename and
     * will create `--version.*.telomere.bed` side-effects in the workdir.
     */
    tuple val("${task.process}"), val('find_telomere'), eval("find_telomere 2>&1 | head -n 1 || true"), topic: versions, emit: versions_findtelomere

    when:
    task.ext.when == null || task.ext.when

    script:

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FINDTELOMERE module does not support Conda. Please use Docker / Singularity instead."
    }

    def prefix          = task.ext.prefix ?: "${meta.id}"
    def window_args     = task.ext.args ?: '99.9 0.1'
    def split_opt       = split_windows ? '--split ' : ''
    def stdout_redirect = split_windows ? '' : "> ${prefix}.windows"
    def max_heap_size_mega = (task.memory.toMega() * 0.9).intValue()
    def max_stack_size_mega = 999
    """
    find_telomere $reference $telomereseq | awk '{print \$1"\\t"\$(NF-4)"\\t"\$(NF-3)"\\t"\$(NF-2)"\\t"\$(NF-1)"\\t"\$NF}' - > ${prefix}.telomere

    java \\
        -Xmx${max_heap_size_mega}M \\
        -Xss${max_stack_size_mega}M \\
        -cp "\$TELOMERE_JAR" \\
        FindTelomereWindows \\
        ${split_opt}${prefix}.telomere \\
        ${window_args} \\
        ${stdout_redirect}
    """

    stub:

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FINDTELOMERE module does not support Conda. Please use Docker / Singularity instead."
    }

    def prefix          = task.ext.prefix ?: "${meta.id}"
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
