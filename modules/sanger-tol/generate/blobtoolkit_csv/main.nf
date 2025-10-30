process GENERATE_BLOBTOOLKIT_CSV {
    tag "${meta.id}"
    label "process_low"
    executor 'local'

    input:
    tuple val(meta), path(pacbio_files)
    val(library_layout)

    output:
    tuple val(meta),    path("samplesheet.csv") , emit: csv
    path("versions.yml")                        , emit: versions

    exec:
    // Note: Manually bump version number when updating module
    def VERSION = "1.0.0"

    def fasta_files = pacbio_files.findAll { it.name.endsWith('.fasta.gz') || it.name.endsWith('.fa.gz') }

    assert fasta_files.size() > 0
    def samplesheet_entries = ["sample,datatype,datafile,library_layout"]
    fasta_files.withIndex().each{ fa, idx -> samplesheet_entries << "${meta.id}_T${idx+1},pacbio,${fa},${library_layout}" }
    file(task.workDir.resolve("samplesheet.csv")).text = samplesheet_entries.join("\n")

    file("${task.workDir}/versions.yml").text = """\
        GENERATE_BLOBTOOLKIT_CSV:
            generate_blobtoolkit_csv: ${VERSION}
        """.stripIndent()

}
