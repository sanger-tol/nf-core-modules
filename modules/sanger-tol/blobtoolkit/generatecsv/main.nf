process BLOBTOOLKIT_GENERATECSV {
    tag "${[meta, meta2, meta3].collect { it.id }.findAll().unique()[0]}"

    input:
    // It is expected that each meta will contain
    // meta.library_layout : {true/false}
    // meta.datatype: {pacbio/ont/illumina}
    tuple val(meta),    path(pacbio_files)
    tuple val(meta2),   path(ont_files)
    tuple val(meta3),   path(illumina_files)

    output:
    tuple val(new_meta),    path("${prefix}.samplesheet.csv")   , emit: csv
    path("versions.yml")                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    exec:
    // Note: Manually bump version number when updating module
    def VERSION       = "1.0.0"

    // CHECK META.ID, IF 1 NOT SET THEN TRY NEXT
    def sample_id = meta.id ?: meta2.id ?: meta3.id
    prefix = task.ext.prefix ?: sample_id

    // SET NEW META OBJECT
    new_meta = meta ?: meta2 ?: meta3

    // CHECK ALL FOR DATA
    def pacbio_data   = pacbio_files.findAll { it.name.endsWith('.fasta.gz') || it.name.endsWith('.fa.gz') }
    def ont_data      = ont_files.findAll { it.name.endsWith('.fasta.gz') || it.name.endsWith('.fa.gz') }
    def illumina_data = illumina_files.findAll { it.name.endsWith('.fasta.gz') || it.name.endsWith('.fa.gz') }

    // ENSURE WE HAVE AT LEAST 1 FILE
    assert (pacbio_data.size() + ont_data.size() + illumina_data.size()) > 0

    // SET HEADER
    def samplesheet_entries = [["sample", "datatype", "datafile", "library_layout"].join(",")]

    // START TOTAL FILE COUNTER
    def counter = 0

    // FORCE PREFIX FOR ALL THE BELOW TO ENSURE SAME NAME
    pacbio_data.eachWithIndex { fa, idx ->
        samplesheet_entries << "${prefix}_T${idx+1+counter},${meta.data_type},${fa},${meta.library_layout}"
    }
    counter += pacbio_data.size()

    ont_data.eachWithIndex { fa, idx ->
        samplesheet_entries << "${prefix}_T${idx+1+counter},${meta2.data_type},${fa},${meta2.library_layout}"
    }
    counter += ont_data.size()

    illumina_data.eachWithIndex { fa, idx ->
        samplesheet_entries << "${prefix}_T${idx+1+counter},${meta3.data_type},${fa},${meta3.library_layout}"
    }

    file(task.workDir.resolve("${prefix}.samplesheet.csv")).text = samplesheet_entries.join("\n")

    file("${task.workDir}/versions.yml").text = """\
        BLOBTOOLKIT_GENERATECSV:
            blobtoolkit_generatecsv: ${VERSION}
        """.stripIndent()
}
