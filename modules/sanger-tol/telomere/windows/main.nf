process TELOMERE_WINDOWS {
    tag "${meta.id}"
    label 'process_low'

    container 'sanger-tol/telomere:0.0.1-c3'

    input:
    tuple val(meta), path(telomere)
    val split

    output:
    tuple val(meta), path("*.windows") , emit: windows
    tuple val("${task.process}"), val('find_telomere_windows'), val("1.0.0"), topic: versions, emit: versions_telomerewindows

    when:
    task.ext.when == null || task.ext.when

    script:
   
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def jar_override = task.ext.telomere_jar ?: ''
    def thresholds  = task.ext.args ?: '99.9 0.1'
    def split_opt   = split ? '--split ' : ''
    // Non-split mode prints windows to stdout; --split writes *.fwd.windows / *.rev.windows (see FindTelomereWindows.java)
    def stdout_redirect = split ? '' : "> ${prefix}.windows"

    // Dynamically generate java mem needs based on task.memory
    // Taken from: nf-core/umicollapse
    def max_heap_size_mega = (task.memory.toMega() * 0.9).intValue()
    def max_stack_size_mega = 999 //most java jdks will not allow Xss > 1GB, so fixing this to the allowed max

    """
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
        echo "TELOMERE_WINDOWS: could not locate telomere.jar (set process ext.telomere_jar or VGP_PIPELINE)" >&2
        exit 1
    fi

    java \\
        -Xmx${max_heap_size_mega}M \\
        -Xss${max_stack_size_mega}M \\
        -cp "\$TELOMERE_JAR" \\
        FindTelomereWindows \\
        ${split_opt}${telomere} \\
        ${thresholds} \\
        ${stdout_redirect}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if ${split}; then
        touch ${prefix}.fwd.windows ${prefix}.rev.windows
    else
        touch ${prefix}.windows
    fi
    """

}
