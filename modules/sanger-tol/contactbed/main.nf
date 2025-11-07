process CONTACTBED {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
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
    // WARNING: This module ships bed_to_contacts.sh as a module binary in
    // ${moduleDir}/resources/usr/bin. Ensure either nextflow.enable.moduleBinaries = true
    // or copy the script into ${projectDir}/bin before using this module.
    """
    bed_to_contacts.sh ${file} > ${meta.id}_contacts.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bed_to_contacts: \$(bed_to_contacts.sh -v)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_contacts.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bed_to_contacts: \$(bed_to_contacts.sh -v)
    END_VERSIONS
    """
}
