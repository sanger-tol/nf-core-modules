process CONTACTBED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}_contacts.bed"
    """
    paste -d '\t' - - < ${file} \
      | awk 'BEGIN {FS="\t"; OFS="\t"} {if (\$1 > \$7) {print substr(\$4,1,length(\$4)-2),\$12,\$7,\$8,"16",\$6,\$1,\$2,"8",\$11,\$5} else {print substr(\$4,1,length(\$4)-2),\$6,\$1,\$2,"8",\$12,\$7,\$8,"16",\$5,\$11} }' \
      | tr '\\-+' '01'  \
      | LC_ALL=C sort -k3,3d -k7,7d \
      | awk 'NF==11' \
      > ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paste: "coreutils 9.5"
        awk: "coreutils 9.5"
        sort: "coreutils 9.5"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_contacts.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paste: "coreutils 9.5"
        awk: "coreutils 9.5"
        sort: "coreutils 9.5"
    END_VERSIONS
    """
}
