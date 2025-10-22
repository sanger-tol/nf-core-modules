//
// MODULE IMPORT BLOCK
//
include { TELOMERE_REGIONS          } from '../../../modules/sanger-tol/telomere/regions/main'
include { GAWK                      } from '../../../modules/nf-core/gawk/main'
include { TELOMERE_WINDOWS          } from '../../../modules/sanger-tol/telomere/windows/main'
include { TELOMERE_EXTRACT          } from '../../../modules/sanger-tol/telomere/extract/main'

workflow TELO_FINDER {

    take:
    ch_reference        // Channel [ val(meta), path(fasta) ]
    ch_telomereseq      // Channel.of( telomere sequence )
    val_split_telomere  // bool

    main:
    ch_versions         = Channel.empty()


    //
    // MODULE: FINDS THE TELOMERIC SEQEUNCE IN REFERENCE
    //
    TELOMERE_REGIONS (
        ch_reference,
        ch_telomereseq
    )
    ch_versions         = ch_versions.mix(TELOMERE_REGIONS.out.versions)


    //
    // MODULE: SPLIT THE TELOMERE FILE INTO 5' and 3' FILES
    //
    if (val_split_telomere) {

        GAWK (
            TELOMERE_REGIONS.out.telomere,
            [],
            true
        )
        ch_versions     = ch_versions.mix(GAWK.out.versions)

        //
        // LOGIC: COLLECT FILES AND ITERATE THROUGH
        //          CHANGE meta.id BASED ON THE FILE NAME
        //          THIS PRODUCES A TRIO OF CHANNELS: [meta], file
        //          THESE CAN BE 0 SIZE AND CAUSE ISSUES DOWNSTREAM
        //
        GAWK.out.output
            .flatMap { meta, files ->
                files
                    .findAll { file -> file.size() > 0 }
                    .collect { file ->
                        if (file.name.contains("direction.0")) {
                            new_meta = meta + [direction: "5P"]
                        }
                        if (file.name.contains("direction.1")) {
                            new_meta = meta + [direction: "3P"]
                        }
                        [new_meta, file]
                    }
            }
            .mix(TELOMERE_REGIONS.out.telomere)
            .set { ch_regions_for_extraction }


    } else {
        ch_regions_for_extraction  = TELOMERE_REGIONS.out.telomere
    }


    //
    // MODULE: GENERATES A WINDOWS FILE FROM THE ABOVE
    //
    TELOMERE_WINDOWS (
        ch_regions_for_extraction
    )
    ch_versions         = ch_versions.mix(TELOMERE_WINDOWS.out.versions)


    //
    // LOGIC: OUTPUT CAN HAVE SIZE 0 WHICH BREAKS gawk IN EXTRACT
    //          FILTER OUT THE 0 SIZE FILES
    //
    TELOMERE_WINDOWS.out.windows
        .filter { _meta, file ->
            file.size() > 0
        }
        .set { ch_filtered_windows_for_extraction  }

    //
    // MODULE: Extract the telomere data from the FIND_TELOMERE
    //          file and reformat into bed
    //
    TELOMERE_EXTRACT(
        ch_filtered_windows_for_extraction
    )
    ch_versions         = ch_versions.mix(TELOMERE_EXTRACT.out.versions)


    TELOMERE_EXTRACT.out.bedgraph
        .map { meta, bedgraph ->
            [ meta - meta.subMap("direction"), bedgraph ]
        }
        .groupTuple(by: 0)
        .set { ch_telo_bedgraphs }


    emit:
    bedgraph_file       = ch_telo_bedgraphs    // Used in pretext_graph
    versions            = ch_versions
}
