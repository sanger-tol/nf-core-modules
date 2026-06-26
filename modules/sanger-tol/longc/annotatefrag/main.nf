process LONGC_ANNOTATEFRAG {
    tag "${meta.id}"
    label 'process_medium'

    // Note: the versions here need to match the versions used in the Wave container below and environment.yml
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/pysam:0.23.0--py312h47d5410_0' :
        'quay.io/biocontainers/pysam:0.23.0--py312h47d5410_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bam.bai"), emit: index
    tuple val("${task.process}"), val('python'), eval('python -c "import sys; print(sys.version.split()[0])"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('annotate_frag.py'), eval("annotate_frag.py --version | sed 's/^.* //'"), topic: versions, emit: versions_annotate_frag

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARNING: This module includes annotate_frag.py as a module binary in
    // ${moduleDir}/resources/usr/bin/annotate_frag.py. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or set the option
    // nextflow.enable.moduleBinaries = true in your nextflow.config file.
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def min_mapq = task.ext.min_mapq ?: 0
    def threads = task.ext.threads ?: task.cpus
    """
    annotate_frag.py \\
        --input ${bam} \\
        --output ${prefix}_annotated.bam \\
        --min-mapq ${min_mapq} \\
        --threads ${threads} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.bam
    touch ${prefix}_annotated.bam.bai
    """
}
