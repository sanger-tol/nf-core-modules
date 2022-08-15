process PRETEXT_MAP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::mulled-v2-f3591ce8609c7b3b33e5715333200aa5c163aa61:c3fe0c6711d0062519837db1b50d17f6899f65e4-0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f3591ce8609c7b3b33e5715333200aa5c163aa61:c3fe0c6711d0062519837db1b50d17f6899f65e4-0' :
        'quay.io/biocontainers/mulled-v2-f3591ce8609c7b3b33e5715333200aa5c163aa61:c3fe0c6711d0062519837db1b50d17f6899f65e4-0' }"

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
