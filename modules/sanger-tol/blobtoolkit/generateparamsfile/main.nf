process BLOBTOOLKIT_GENERATEPARAMSFILE {
    tag "${meta.id}"
    executor "local"

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
    tuple val(meta), path("*.json") , emit: json_params_file
    tuple val("${task.process}"), val('blobtoolkit_generateparamsfile'), val('1.0.0'), emit: versions_blobtoolkitgenerateparamsfile, topic: versions

    when:
    task.ext.when == null || task.ext.when

    exec:
    prefix = task.ext.prefix ?: "${meta.id}"

    def btk_inputs = btk_extra_opts + [
        'fasta': reference.toUriString(),
        'busco_lineages': busco_lineages,
        'taxon': taxon,
        'taxdump': tax_dump?.toUriString(),
        'blastp': blastp?.toUriString(),
        'blastn': blastn?.toUriString(),
        'blastx': blastx?.toUriString(),
    ].findAll { pair -> pair.value } // filter out falsy values (null, false, "", [], etc)

    def jsonBuilder = new groovy.json.JsonBuilder(btk_inputs)
    file("${task.workDir}/${prefix}.blobtoolkit_params_file.json").text = jsonBuilder.toPrettyString()
}
