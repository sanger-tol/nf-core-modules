include { PRETEXTMAP        } from '../../../modules/nf-core/pretextmap/main'
include { PRETEXTSNAPSHOT   } from '../../../modules/nf-core/pretextsnapshot/main'
// include { PRETEXTANNOTATE   } from '../../../modules/local/pretextannotate/main'

workflow PRETEXTMAP_SNAPSHOT_ANNOTATE {
    take:
    ch_reference    // [[meta], reference, fai  ]
    ch_mapped_bam   // [[meta], mapped_bam      ]
    // ch_chromlist

    main:
    ch_versions = channel.empty()

    //
    // MODULE: PROCESS MAPPED BAMS INTO PRETEXTMAP
    //
    PRETEXTMAP (
        ch_mapped_bam,
        ch_reference
    )


    //
    // MODULE: GRAB PNG OF INPUT PRETEXT FILE
    //
    PRETEXTSNAPSHOT (
        PRETEXTMAP.out.pretext
    )


    //
    // MODULE: ANNOTATE THE PNG WITH CHROMOSOME DATA
    //         This is to be added in the future!
    //         The subworkflow needs adding in before I have
    //         time to write the script.
    //
    // PRETEXTANNOTATE (
    //      PRETEXTSNAPSHOT.out.image,
    //      ch_chromlist
    // )


    emit:
    pretext_file = PRETEXTMAP.out.pretext
    pretext_pngs = PRETEXTSNAPSHOT.out.image
    // pretext_annoated_image = PRETEXTANNOTATE.out.png
}
