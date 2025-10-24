process SAMTOOLS_GREPREADGROUPS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0':
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), stdout, emit: read_group
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -H ${bam} |\\
        ## Do some pretty printing with awk to avoid starting or trailing newlines
        awk '/@RG/ { printf "%s%s", (n++ ? "\\n" : ""), \$0 }'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -n "@RG	ID:49240_6#4	DT:2024-07-25T16:52:42Z	PU:20240725_LH00275_0063_A22CHVNLT4_6#4	LB:72611119	PG:SCS	SM:SAMEA111361848	CN:SC	PL:ILLUMINA	DS:ERP129860: Sequencing and assembly of symbiont/host genomes in marine and freshwater species, where at least one partner is a microbe, for the Aquatic Symbiosis Genomics project. Upon completion, genomes will be publicly available on an online data coordination center run by the European Bioinformatics Institute. All the data the project generates will be released openly and fairly, as the effort works towards a global community of scientists and ecologists aiming to understand and conserve ocean and freshwater biodiversity. This data is part of a pre-publication release. For information on the proper use of pre-publication data shared by the Wellcome Trust Sanger Institute (including details of any publication moratoria), please see http://www.sanger.ac.uk/datasharing/"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
