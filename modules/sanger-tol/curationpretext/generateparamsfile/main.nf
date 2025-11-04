process CURATIONPRETEXT_GENERATEPARAMSFILE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reference)
    val longread_dir
    val cram_dir
    val telomere_motif
    val aligner
    val cpretext_extra_opts

    output:
    tuple val(meta), path("${prefix}.curationpretext_params_file.json"), emit: params_file
    path("versions.yml")                                               , emit: versions

    exec:
    def VERSION = "1.0.0"

    prefix = task.ext.prefix ?: "${meta.id}"

    def cpretext_inputs = cpretext_extra_opts + [
        'sample': meta.id,
        'reads': longread_dir.toUriString(),
        'cram': cram_dir.toUriString(),
        'teloseq': telomere_motif,
        'aligner': aligner,
    ].findAll { it.value } // filter out falsy values (null, false, "", [], etc)
    def jsonBuilder = new groovy.json.JsonBuilder(cpretext_inputs)
    file("${task.workDir}/cpretext_params_file.json").text = jsonBuilder.toPrettyString()

    file("${task.workDir}/versions.yml").text = """\
        CURATIONPRETEXT_GENERATEPARAMSFILE:
            curationpretext_generateparamsfile: ${VERSION}
        """.stripIndent()
}
