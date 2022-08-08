#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { SELFCOMP_SPLITFASTA } from '../../../../modules/selfcomp/splitfasta/main.nf'
workflow test_selfcomp_splitfasta {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.tol_test_data['small_genome']['pEimTen1']['assembly']['canu_contigs_fasta'], checkIfExists: true)
    ]
    SELFCOMP_SPLITFASTA ( input )
}