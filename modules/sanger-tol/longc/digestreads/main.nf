process LONGC_DIGESTREADS {
    tag "${meta.id}"
    label 'process_medium'

    // Note: the versions here need to match the versions used in the Wave container below and environment.yml
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a7/a7968a0e334fc598eebd9dede58fff52ca342b26697241222a7681171447c992/data' :
        'community.wave.seqera.io/library/samtools_python:40e23973bbc3d3dd' }"

    input:
    tuple val(meta), path(reads), path(reference)
    val cutter    // restriction enzyme name, e.g. 'NlaIII'

    output:
    tuple val(meta), path("*.fastq.gz"), emit: digested_reads
    tuple val("${task.process}"), val('python'), eval('python -c "import sys; print(sys.version.split()[0])"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('digest_reads.py'), eval("digest_reads.py --version | sed 's/^.* //'"), topic: versions, emit: versions_digest_reads
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n1 | sed "s/ /,/"'), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARNING: This module includes digest_reads.py as a module binary in
    // ${moduleDir}/resources/usr/bin/digest_reads.py. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or set the option
    // nextflow.enable.moduleBinaries = true in your nextflow.config file.
    //
    // reference is required for CRAM; pass assembly FASTA for all input types.
    def prefix   = task.ext.prefix   ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def reads_list = reads instanceof List ? reads : [ reads ]
    def reads_arg = reads_list.join(' ')
    def is_cram = reads_list.every { n -> n.name.toLowerCase().endsWith('.cram') }
    def is_bam = reads_list.every { n -> n.name.toLowerCase().endsWith('.bam') }
    def is_alignment = is_cram || is_bam
    def digest_in = is_alignment ? '-' : reads_arg
    def samtools_ref = is_cram ? "--reference ${reference}" : ''
    def samtools_pipe = is_alignment ? "samtools fastq -T '*' --threads ${task.cpus} ${samtools_ref} ${reads_arg} |" : ''
    """
    ${samtools_pipe} \\
    digest_reads.py \\
        --cutter ${cutter} \\
        ${args} \\
        ${digest_in} | \\
    python -c "import gzip, sys; gzip.open('${prefix}.fastq.gz', 'wt').write(sys.stdin.read())"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -c "import gzip; gzip.open('${prefix}.fastq.gz', 'wb').write(b'')"
    """
}
