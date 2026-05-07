process FINDTELOMERE {
    tag "${meta.id}"
    label 'process_low'

    container 'sanger-tol/telomere:0.0.2-c1'

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

    when:
    task.ext.when == null || task.ext.when

    script:

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FINDTELOMERE module does not support Conda. Please use Docker / Singularity instead."
    }

    def prefix          = task.ext.prefix ?: "${meta.id}"
    def jar_override    = task.ext.telomere_jar ?: ''
    def window_args     = task.ext.args ?: '99.9 0.1'
    def split_opt       = split_windows ? '--split ' : ''
    def stdout_redirect = split_windows ? '' : "> ${prefix}.windows"
    def max_heap_size_mega = (task.memory.toMega() * 0.9).intValue()
    def max_stack_size_mega = 999
    """
    find_telomere $reference $telomereseq | awk '{print \$1"\\t"\$(NF-4)"\\t"\$(NF-3)"\\t"\$(NF-2)"\\t"\$(NF-1)"\\t"\$NF}' - > ${prefix}.telomere

    TELOMERE_JAR="${jar_override}"
    if [ -z "\$TELOMERE_JAR" ]; then
        if [ -f telomere.jar ]; then
            TELOMERE_JAR=telomere.jar
        elif [ -f /opt/telomere/telomere.jar ]; then
            TELOMERE_JAR=/opt/telomere/telomere.jar
        else
            _ft=\$(command -v find_telomere 2>/dev/null || true)
            if [ -n "\$_ft" ]; then
                _bin=\$(readlink -f "\$_ft" 2>/dev/null || readlink "\$_ft" 2>/dev/null || echo "\$_ft")
                _dir=\$(dirname "\$_bin")
                if [ -f "\$_dir/telomere.jar" ]; then
                    TELOMERE_JAR="\$_dir/telomere.jar"
                fi
            fi
        fi
        if [ -z "\$TELOMERE_JAR" ] && [ -n "\${VGP_PIPELINE:-}" ] && [ -f "\${VGP_PIPELINE}/telomere/telomere.jar" ]; then
            TELOMERE_JAR="\${VGP_PIPELINE}/telomere/telomere.jar"
        fi
        if [ -z "\$TELOMERE_JAR" ] && [ -f /usr/local/share/telomere/telomere.jar ]; then
            TELOMERE_JAR=/usr/local/share/telomere/telomere.jar
        fi
    fi
    if [ -z "\$TELOMERE_JAR" ] || [ ! -f "\$TELOMERE_JAR" ]; then
        echo "FINDTELOMERE: could not locate telomere.jar for FindTelomereWindows (set process ext.telomere_jar or VGP_PIPELINE)" >&2
        exit 1
    fi

    # java on PATH comes from the container image; failures are usually wrong/missing -cp (telomere.jar), not missing JVM
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
    def bed_stem        = reference.getName()
    """
    touch ${prefix}.telomere
    touch ${bed_stem}.fwd.telomere.bed ${bed_stem}.rev.telomere.bed
    if ${split_windows}; then
        touch ${prefix}.fwd.windows ${prefix}.rev.windows
    else
        touch ${prefix}.windows
    fi
    """

}
