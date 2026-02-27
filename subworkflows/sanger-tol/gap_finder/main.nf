//
// GENERATE BED FILE OF GAPS AND LENGTH IN REFERENCE
//

include { SEQTK_CUTN        } from '../../../modules/nf-core/seqtk/cutn/main'
include { GAWK              } from '../../../modules/nf-core/gawk/main'
include { TABIX_BGZIPTABIX  } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow GAP_FINDER {
    take:
    ch_reference     // Channel [ val(meta), path(fasta) ]
    val_run_bgzip    // val(boolean)

    main:

    //
    // MODULE: GENERATES A GAP SUMMARY FILE
    //
    SEQTK_CUTN (
        ch_reference
    )


    // NOTE: BUILD THE GAWK FILE
    ch_reorder_gap_file_awk = channel.of('''\
    BEGIN { OFS = "\\t" } {
        print $0, sqrt(($3-$2)*($3-$2))
    }'''.stripIndent())
        .collectFile(name: "reorder_gap_file_awk.awk", cache: true)
        .collect()

    //
    // MODULE: ADD THE LENGTH OF GAP TO BED FILE - INPUT FOR PRETEXT MODULE
    //
    GAWK (
        SEQTK_CUTN.out.bed,
        ch_reorder_gap_file_awk,
        false
    )


    //
    // MODULE: BGZIP AND TABIX THE GAP FILE
    //
    TABIX_BGZIPTABIX (
        SEQTK_CUTN.out.bed.filter{ meta, file -> val_run_bgzip}
    )

    emit:
    gap_file        = GAWK.out.output
    gap_tabix       = TABIX_BGZIPTABIX.out.gz_index
}
