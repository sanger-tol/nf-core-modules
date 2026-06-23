#!/usr/bin/env nextflow

//
// SANGER_TOL SUBWORKFLOW IMPORT BLOCK
//
include { GAP_FINDER                        } from '../../sanger-tol/gap_finder/main'
include { TELO_FINDER                       } from '../../sanger-tol/telo_finder/main'
include { READ_COVERAGE                     } from '../../sanger-tol/read_coverage/main'
include { REPEAT_DENSITY                    } from '../../sanger-tol/repeat_density/main'

workflow ACCESSORY_FILES {
    take:
    ch_reference_tuple      // Channel [ val(meta), path(file)   ]
    ch_reference_sizes      // Channel [ val(meta), path(file)   ]
    ch_longread_reads       // Channel [ val(meta), [path(file)] ]
    val_teloseq             // val(telomere_sequence)
    val_run_gap_finder      // val(bool)
    val_run_telo_finder     // val(bool)
    val_run_repeat_density  // val(bool)
    val_run_coverage        // val(bool)
    val_split_telomere      // val(bool)


    main:

    //
    // SUBWORKFLOW: GENERATES A GAP.BED FILE TO ID THE LOCATIONS OF GAPS
    //
    GAP_FINDER (
        reference_tuple.filter{ _meta, _file -> run_gap_finder },
        false
    )


    //
    // SUBWORKFLOW: GENERATE TELOMERE WINDOW FILES WITH LONGREAD READS AND REFERENCE
    //
    TELO_FINDER (
        reference_tuple.filter{ _meta, _file -> run_telo_finder },
        val_teloseq,
        val_split_telomere,
        false
    )


    //
    // SUBWORKFLOW: GENERATES A BIGWIG FOR A REPEAT DENSITY TRACK
    //
    REPEAT_DENSITY (
        reference_tuple.filter{ _meta, _file -> run_repeat_density },
        ch_reference_sizes
    )


    //
    // SUBWORKFLOW: Takes reference, longread reads
    //
    READ_COVERAGE (
        longread_reads,
        reference_tuple.filter{ _meta, _file -> run_coverage },
        ch_reference_sizes.map{ _meta, file -> file }
    )


    emit:
    gap_file        = GAP_FINDER.out.gap_file
    repeat_file     = REPEAT_DENSITY.out.repeat_density
    telo_file       = TELO_FINDER.out.bedgraph_file      // This is the possible collection of telomere files
    coverage_output = READ_COVERAGE.out.bigwig
}
