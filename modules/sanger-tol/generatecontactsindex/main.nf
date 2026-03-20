process GENERATE_CONTACTS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(contacts)

    output:
    tuple val(meta), path(contacts), path("*.index.tsv"), emit: contacts_with_index
    tuple val("${task.process}"), val('coreutils'), eval('ls --version | sed -n "s/ls (GNU coreutils) //p"'), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create index TSV file from contacts
    # Index typically contains unique contact pairs or sorted positions
    cut -f1,2 ${contacts} |\\
    LC_ALL=C sort \\
        -k1,1 -k2,2n \\
        -u -S ${task.memory.toGiga()}G \\
        > ${prefix}.index.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.index.tsv
    """
}
