//
// MODULE IMPORT BLOCK
//
include { TELOMERE_REGIONS          } from '../../../modules/sanger-tol/telomere/regions/main'
//include { GAWK_SPLIT_DIRECTIONS     } from '../../../modules/local/gawk_split_directions/main'
include { GAWK} from '../../../modules/nf-core/gawk/main'
include { TELOMERE_WINDOWS          } from '../../../modules/sanger-tol/telomere/windows/main'
include { TELOMERE_EXTRACT          } from '../../../modules/sanger-tol/telomere/extract/main'

workflow TELO_FINDER {

    take:
    ch_reference     // Channel [ val(meta), path(fasta) ]
    telomereseq      // Channel.of( telomere sequence )
    split_telomere   // bool

    main:
    ch_versions         = Channel.empty()


    //
    // MODULE: FINDS THE TELOMERIC SEQEUNCE IN REFERENCE
    //
    TELOMERE_REGIONS (
        ch_reference,
        telomereseq
    )
    ch_versions         = ch_versions.mix(TELOMERE_REGIONS.out.versions)


    //
    // MODULE: SPLIT THE TELOMERE FILE INTO 5' and 3' FILES
    //              THIS IS RUNNING ON A LOCAL VERSION OF THE GAWK MODULE
    //
    if (split_telomere) {

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
                        def new_meta = meta.clone()
                        if (file.name.contains("direction.0")) {
                            new_meta.id = "${meta.id}_5P"
                        }
                        if (file.name.contains("direction.1")) {
                            new_meta.id = "${meta.id}_3P"
                        }
                        [new_meta, file]
                    }
            }
            .mix(TELOMERE_REGIONS.out.telomere)
            .set { for_extraction }


    } else {
        for_extraction  = TELOMERE_REGIONS.out.telomere
    }


    //
    // MODULE: GENERATES A WINDOWS FILE FROM THE ABOVE
    //
    TELOMERE_WINDOWS (
        for_extraction
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
        .set { safe_extract_input  }

    //
    // MODULE: Extract the telomere data from the FIND_TELOMERE
    //          file and reformat into bed
    //
    TELOMERE_EXTRACT(
        safe_extract_input
    )
    ch_versions         = ch_versions.mix(TELOMERE_EXTRACT.out.versions)


    TELOMERE_EXTRACT.out.bedgraph
        .map{ _meta, bedgraph ->
            bedgraph
        }
        .collect()
        .set { telo_bedgraphs }


    emit:
    bedgraph_file       = telo_bedgraphs    // Used in pretext_graph
    versions            = ch_versions
}
