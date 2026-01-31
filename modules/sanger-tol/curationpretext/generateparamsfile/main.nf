process CURATIONPRETEXT_GENERATEPARAMSFILE {
    tag "${meta.id}"
    executor "local"

    input:
    tuple val(meta), path(reference)
    val longread_dir
    val cram_dir
    val telomere_motif
    val aligner
    val cpretext_extra_opts

    output:
    tuple val(meta), path("${prefix}.curationpretext_params_file.json"), emit: json_params_file
    path("versions.yml")                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    exec:
    def VERSION = "1.0.0"

    prefix = task.ext.prefix ?: "${meta.id}"

    def cpretext_inputs = cpretext_extra_opts + [
        'sample': meta.id,
        'reads': longread_dir.toUriString(),
        'cram': cram_dir.toUriString(),
        'teloseq': telomere_motif,
        'aligner': aligner,
    ].findAll { kv -> kv.value } // filter out falsy values (null, false, "", [], etc)
    def jsonBuilder = new groovy.json.JsonBuilder(cpretext_inputs)
    file("${task.workDir}/${prefix}.curationpretext_params_file.json").text = jsonBuilder.toPrettyString()

    file("${task.workDir}/versions.yml").text = """\
        CURATIONPRETEXT_GENERATEPARAMSFILE:
            curationpretext_generateparamsfile: ${VERSION}
        """.stripIndent()
}
