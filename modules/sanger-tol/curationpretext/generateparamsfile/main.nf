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
    tuple val(meta), path("*.json"), emit: json_params_file
    tuple val("${task.process}"), val('curationpretext_generateparamsfile'), val('1.0.0'), emit: versions_curationpretextgenerateparamsfile, topic: versions

    when:
    task.ext.when == null || task.ext.when

    exec:
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
}
