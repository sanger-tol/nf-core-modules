process SANGERTOL_GENERATEBLOBTOOLKITCSV {
    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path(input_files)
    val(longread_technology)
    val(library_layout)

    output:
    tuple val(meta),    path("${prefix}.samplesheet.csv")   , emit: csv
    path("versions.yml")                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    exec:
    // Note: Manually bump version number when updating module
    def VERSION     = "1.0.0"
    prefix          = task.ext.prefix ?: "${meta.id}"

    def fasta_files = input_files.findAll { it.name.endsWith('.fasta.gz') || it.name.endsWith('.fa.gz') }

    assert fasta_files.size() > 0
    def samplesheet_entries = ["sample,datatype,datafile,library_layout"]
    fasta_files.withIndex().each{ fa, idx -> samplesheet_entries << "${meta.id}_T${idx+1},${longread_technology},${fa},${library_layout}" }
    file(task.workDir.resolve("${prefix}.samplesheet.csv")).text = samplesheet_entries.join("\n")

    file("${task.workDir}/versions.yml").text = """\
        SANGERTOL_GENERATEBLOBTOOLKITCSV:
            sangertol_generateblobtoolkitcsv: ${VERSION}
        """.stripIndent()

}
