process BUSCOANNOTATION_ASSIGNANCESTRAL {
    tag "$meta.id"
    label 'process_medium'

    def version = '0.001-c3'
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the buscoannotation process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/buscogeneannotation:${version}"

    input:
    tuple val(meta), path(paintedtsv)
    path fulltable



    output:
    tuple val(meta), path("*.csv")                , emit: ancestralresult
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    assign_anc.py \\
        $args \\
        ${paintedtsv} \\
        ${fulltable} \\
        ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        buscoannotation: $version
    END_VERSIONS
    """
}

