#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCOANNOTATION_BUSCOPAINTERLEP } from '../../../../modules/buscoannotation/buscopainterlep/main.nf'

workflow test_buscoannotation_buscopainterlep {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.tol_test_data['small_genome']['Leptidea_sinapis']['genomic_data']['ilLepSina1']['busco']['full_table'], checkIfExists: true)
    ]

    merian = file(params.tol_test_data['small_genome']['Leptidea_sinapis']['genomic_data']['ilLepSina1']['busco']['ancestral_units'], checkIfExists: true) 

    BUSCOANNOTATION_BUSCOPAINTERLEP ( input, merian )
}
