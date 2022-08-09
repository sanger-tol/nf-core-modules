process PRETEXT_MAP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pretextmap:0.1.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pretextmap:0.1.9--h9f5acd7_1' :
        'quay.io/biocontainers/pretextmap:0.1.9--h9f5acd7_1' }"

    input:
    tuple val(meta), path(alignment)
    path(fai)

    output:
    tuple val(meta), path("*pretext"), emit: pretext
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (alignment.split({'.'})[-1] in ['sam','bam']) {
        """
        samtools view -h | PretextMap $args -o $meta.id".pretext";

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            PretextMap: \$(PretextMap | grep 'Version' | cut -f3 -d' ')
            samtools: \$(samtools 2>&1 | grep Version | cut -f2 -d' ')
        END_VERSIONS
        """
    } else {
        // If not a sam/bam file assume it's from juicer pre
        """
        prepare_pretext.sh $fai $alignment | PretextMap $args -o $meta.id".pretext";
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            PretextMap: \$(PretextMap | grep 'Version' | cut -f3 -d' ')
        END_VERSIONS
        """
    }
}
