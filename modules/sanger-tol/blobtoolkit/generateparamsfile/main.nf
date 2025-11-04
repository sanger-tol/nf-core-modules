process BLOBTOOLKIT_GENERATEPARAMSFILE {
    tag "${meta.id}"

    input:
    tuple val(meta), val(reference)
    val blastp
    val blastn
    val blastx
    val tax_dump
    val busco_lineages
    val taxon
    val btk_extra_opts

    output:
    tuple val(meta), path("${prefix}.blobtoolkit_params_file.json") , emit: json_params_file
    path("versions.yml")                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    exec:
    // Note: Manually bump version number when updating module
    def VERSION = "1.0.0"

    prefix = task.ext.prefix ?: "${meta.id}"

    def btk_inputs = btk_extra_opts + [
        'fasta': reference.toUriString(),
        'busco_lineages': busco_lineages,
        'taxon': taxon,
        'taxdump': tax_dump?.toUriString(),
        'blastp': blastp?.toUriString(),
        'blastn': blastn?.toUriString(),
        'blastx': blastx?.toUriString(),
    ].findAll { it.value } // filter out falsy values (null, false, "", [], etc)
    def jsonBuilder = new groovy.json.JsonBuilder(btk_inputs)
    file("${task.workDir}/${prefix}.blobtoolkit_params_file.json").text = jsonBuilder.toPrettyString()

    file("${task.workDir}/versions.yml").text = """\
        BLOBTOOLKIT_GENERATEPARAMSFILE:
            blobtoolkit_generateparamsfile: ${VERSION}
        """.stripIndent()
}
