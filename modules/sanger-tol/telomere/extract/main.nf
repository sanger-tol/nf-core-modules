process TELOMERE_EXTRACT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(telomere)

    output:
    tuple val(meta), path("*.bed")      , emit: bed
    tuple val(meta), path("*.bedgraph") , emit: bedgraph
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    awk 'BEGIN { OFS = "\t" }
    {
        gsub(">", "")
        print \$2, \$4, \$5 >> "${prefix}_telomere.bed"
        print \$2, \$4, \$5, (((\$5-\$4)<0)?-(\$5-\$4):(\$5-\$4)) >> "${prefix}_telomere.bedgraph"
    }' ${telomere}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_telomere.bed
    touch ${prefix}_telomere.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
