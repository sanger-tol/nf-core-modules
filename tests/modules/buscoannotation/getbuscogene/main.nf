#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCOANNOTATION_GETBUSCOGENE } from '../../../../modules/buscoannotation/getbuscogene/main.nf'

workflow test_buscoannotation_getbuscogene {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.tol_test_data['small_genome']['Leptidea_sinapis']['genomic_data']['ilLepSina1']['busco']['full_table'], checkIfExists: true)
    ]

    BUSCOANNOTATION_GETBUSCOGENE ( input )
}
