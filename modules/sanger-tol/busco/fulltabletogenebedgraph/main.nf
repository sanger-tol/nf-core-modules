process BUSCO_FULLTABLE_TO_GENE_BEDGRAPH {
    tag "${meta.id}"
    label 'process_single'
     
    conda "${moduleDir}/environment.yml"
    container "docker://amancevice/pandas:latest"

    input:
    tuple val(meta), path(full_table)

    output:
    tuple val(meta), path("complete_buscos.bedgraph"), emit: complete_bedgraph
    tuple val(meta), path("duplicated_buscos.bedgraph"), emit: duplicated_bedgraph
    tuple val(meta), path("fragmented_buscos.bedgraph"), emit: fragmented_bedgraph
    tuple val("${task.process}"), val('python'), eval('python -c "import sys; print(sys.version.split()[0])"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas as pd; print(pd.__version__)"'), emit: versions_pandas, topic: versions
    tuple val("${task.process}"), val('busco_fulltable_to_bedgraph'), eval("python ${moduleDir}/resources/usr/bin/busco_fulltable_to_bedgraph --version | sed 's/^.* //'"), emit: versions_busco_fulltable_to_bedgraph, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    busco_fulltable_to_bedgraph \\
        ${full_table} \\
        -o .
    """

    stub:
    """
    touch complete_buscos.bedgraph
    touch duplicated_buscos.bedgraph
    touch fragmented_buscos.bedgraph
    """
}
