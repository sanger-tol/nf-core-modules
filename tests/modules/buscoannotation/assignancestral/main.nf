#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCOANNOTATION_ASSIGNANCESTRAL } from '../../../../modules/buscoannotation/assignancestral/main.nf'

workflow test_buscoannotation_assignancestral {
    
    ancestral_locations = [
        [ id:'test', single_end:false ],
        file(params.tol_test_data['small_genome']['Oscheius_sp']['genomic_data']['ilLepSina1']['busco']['buscopaintedfile'], checkIfExists: true)
    ]

    fulltable = file(params.tol_test_data['small_genome']['Oscheius_sp']['genomic_data']['ilLepSina1']['busco']['full_table'], checkIfExists: true)

    BUSCOANNOTATION_ASSIGNANCESTRAL ( ancestral_locations, fulltable  )
}
