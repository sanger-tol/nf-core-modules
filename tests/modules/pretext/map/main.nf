#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRETEXT_MAP } from '../../../../modules/pretext/map/main.nf'

workflow test_pretext_map {
    
    input = [
        [ id:'test' ], // meta map
        file(params.tol_test_data['small_genome']['pEimTen1']['assembly']['alignments_sorted'], checkIfExists: true)]

    fai = file(params.tol_test_data['small_genome']['pEimTen1']['assembly']['scaffolds_fai'], checkIfExists: true)
    
    PRETEXT_MAP ( input , fai)
}
