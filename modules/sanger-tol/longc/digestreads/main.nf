process LONGC_DIGESTREADS {
    tag "${meta.id}"
    label 'process_medium'

    // Note: the versions here need to match the versions used in the Wave container below and environment.yml
    // environment.yml / Wave freeze: samtools=1.23.1 + gzip=1.14 + python=3.14.6 + biopython=1.86
    // Docker:  community.wave.seqera.io/library/samtools_gzip_python_biopython:c338e4d96f9312d5
    // Singularity native SIF: samtools_gzip_python_biopython:41a0e9b031779042
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/3867ff7ef51d4be6e84565f1f506f2073ab4027ddd66b79744a534440fc84250/data' :
        'community.wave.seqera.io/library/samtools_gzip_biopython_python:09f17801b9745830' }"

    input:
    tuple val(meta), path(reads, stageAs: 'reads/?/*'), path(reference)
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
    def is_fastq = reads_list.every { n -> n.name.toLowerCase() ==~ /.*\.(fq|fastq|fa|fasta)(\.gz)?$/ }
    def is_alignment = is_cram || is_bam

    if (reads_list.size() > 1 && !(is_alignment || is_fastq)) {
        error("LONGC_DIGESTREADS: all reads in a list must be the same type (FASTQ/FASTA, BAM, or CRAM) for sample ${meta.id}")
    }
    if (is_bam && is_cram) {
        error("LONGC_DIGESTREADS: cannot mix BAM and CRAM files in the same reads list for sample ${meta.id}")
    }

    def samtools_ref = is_cram ? "--reference ${reference}" : ''
    def samtools_pipe = is_alignment ? "samtools fastq -T '*' --threads ${task.cpus} ${samtools_ref} ${reads_arg} |" : ''
    def fastq_cat_pipe = ''
    if (!is_alignment && reads_list.size() > 1) {
        def cat_cmds = reads_list.collect { f ->
            f.name.toLowerCase().endsWith('.gz') ? "gzip -dc '${f}'" : "cat '${f}'"
        }.join('\n    ')
        fastq_cat_pipe = "{\n    ${cat_cmds}\n} |"
    }
    def digest_in = (is_alignment || reads_list.size() > 1) ? '-' : reads_list[0]
    """
    set -o pipefail
    ${samtools_pipe}${fastq_cat_pipe}
    digest_reads.py \\
        --cutter ${cutter} \\
        ${args} \\
        ${digest_in} | \\
    gzip -c > ${prefix}_${cutter}.fastq.gz

    # Fail loudly if digestion produced no monomers (gzip of empty stdin is still a non-zero .gz).
    python3 -c 'import gzip,sys; sys.exit(0 if gzip.open(sys.argv[1],"rb").read(1) else 1)' ${prefix}_${cutter}.fastq.gz \\
        || { echo "ERROR: ${prefix}_${cutter}.fastq.gz is empty after digestion" >&2; exit 1; }
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_${cutter}.fastq.gz
    """
}
