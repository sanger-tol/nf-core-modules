process GENERATE_CONTACTS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(contacts)

    output:
    tuple val(meta), path(contacts), path("*.index.tsv"), emit: contacts_with_index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Create index TSV file from contacts
    # Index typically contains unique contact pairs or sorted positions
    LC_ALL=C sort -k1,1 -k2,2n -S 2G ${contacts} | cut -f1,2 | uniq > ${meta.id}.index.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: "coreutils 9.1"
        cut: "coreutils 9.1"
        uniq: "coreutils 9.1"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.index.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: "coreutils 9.1"
        cut: "coreutils 9.1"
        uniq: "coreutils 9.1"
    END_VERSIONS
    """
}

