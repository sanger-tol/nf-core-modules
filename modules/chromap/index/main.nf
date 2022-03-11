process CHROMAP_INDEX {
    tag '$fasta'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::chromap=0.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chromap:0.2.0--hd03093a_1' :
        'quay.io/biocontainers/chromap:0.2.0--hd03093a_1' }"

    input:
    path fasta

    output:
    path "*.index"     , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = fasta.baseName
    """
    chromap \\
        -i \\
        $args \\
        -t $task.cpus \\
        -r $fasta \\
        -o ${prefix}.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chromap: \$(echo \$(chromap --version 2>&1))
    END_VERSIONS
    """
}
