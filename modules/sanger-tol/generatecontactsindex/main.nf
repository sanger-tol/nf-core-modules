process GENERATE_CONTACTS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e2eac050990ed5d6a558bd4a08c0bf424104db9a1e291e8df24aedf63c0d51b5/data' :
        'community.wave.seqera.io/library/coreutils_gawk:1839a575e817eec0' }"

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
