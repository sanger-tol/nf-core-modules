process BUSCOFULLTABLETOGENEBEDGRAPH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(full_table)

    output:
    tuple val(meta), path("*.complete_buscos.bedgraph"),   emit: complete_bedgraph
    tuple val(meta), path("*.duplicated_buscos.bedgraph"), emit: duplicated_bedgraph
    tuple val(meta), path("*.fragmented_buscos.bedgraph"), emit: fragmented_bedgraph
    tuple val("${task.process}"), val('python'), eval('python -c "import sys; print(sys.version.split()[0])"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas as pd; print(pd.__version__)"'), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('busco_fulltable_to_bedgraph'), eval("python ${moduleDir}/resources/usr/bin/busco_fulltable_to_bedgraph --version | sed 's/^.* //'"), topic: versions, emit: versions_busco_fulltable_to_bedgraph

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    busco_fulltable_to_bedgraph \\
        ${full_table} \\
        -o . \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.complete_buscos.bedgraph
    touch ${prefix}.duplicated_buscos.bedgraph
    touch ${prefix}.fragmented_buscos.bedgraph
    """
}
