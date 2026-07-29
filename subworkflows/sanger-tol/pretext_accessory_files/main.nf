#!/usr/bin/env nextflow

//
// SANGER_TOL SUBWORKFLOW IMPORT BLOCK
//
include { GAP_FINDER                        } from '../gap_finder/main'
include { TELO_FINDER                       } from '../telo_finder/main'
include { READ_COVERAGE                     } from '../read_coverage/main'
include { REPEAT_DENSITY                    } from '../repeat_density/main'

workflow PRETEXT_ACCESSORY_FILES {
    take:
    ch_reference_tuple  // Channel [ val(meta), path(file)   ]
    ch_reference_sizes  // Channel [ val(meta), path(file)   ]
    ch_longread_reads   // Channel [ val(meta), [path(file)] ]
    ch_teloseq          // Channel [ val(meta), path(telomere_sequence) ]
    val_split_telomere  // val(bool)
    val_run_telomere    // val(bool)
    val_run_gaps        // val(bool)
    val_run_coverage    // val(bool)
    val_run_repeat_den  // val(bool)
    val_run_busco_track // val(bool) PLACEHOLDER
    val_run_pebble      // val(bool) PLACEHOLDER


    main:

    //
    // SUBWORKFLOW: GENERATES A GAP.BED FILE TO ID THE LOCATIONS OF GAPS
    //
    GAP_FINDER(
        ch_reference_tuple.filter { val_run_gaps },
        false
    )


    //
    // SUBWORKFLOW: GENERATE TELOMERE WINDOW FILES WITH LONGREAD READS AND REFERENCE
    //
    TELO_FINDER (
        ch_reference_tuple.filter { val_run_telomere },
        ch_teloseq,
        val_split_telomere,
        false
    )

    telomere_data = TELO_FINDER.out.telomere
                        .mix( TELO_FINDER.out.telomere_bed_fwd )
                        .mix( TELO_FINDER.out.telomere_bed_rev )


    //
    // SUBWORKFLOW: GENERATES A BIGWIG FOR A REPEAT DENSITY TRACK
    //
    REPEAT_DENSITY (
        ch_reference_tuple.filter { val_run_repeat_den },
        ch_reference_sizes
    )


    //
    // SUBWORKFLOW: Takes reference, longread reads
    //
    READ_COVERAGE (
        ch_longread_reads,
        ch_reference_tuple.filter { val_run_coverage },
        ch_reference_sizes.map { _meta, file -> file }
    )


    //
    // MODULE: PEBBLE PLACEHOLDER
    //
    // PEBBLE (
    //  ch_reference_tuple
    // )


    //
    // SUBWORKFLOW: BUSCO_TRACK PLACEHOLDER
    //
    // BUSCO_TRACK (
    //  ch_reference_tuple,
    //  busco_data
    //)


    emit:
    gap_file            = GAP_FINDER.out.gap_file
    repeat_file         = REPEAT_DENSITY.out.repeat_density
    telo_file           = telomere_data                     // This is the possible collection of telomere files
    coverage_file       = READ_COVERAGE.out.bigwig
}
