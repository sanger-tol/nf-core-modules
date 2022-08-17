#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SELFCOMP_MUMMER2BED } from '../../../../modules/selfcomp/mummer2bed/main.nf'

workflow test_selfcomp_mummer2bed {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.tol_test_data['small_genome']['Oscheius_sp']['genomic_data']['nxOscSpen2']['selfcomp']['mummer_file'], checkIfExists: true)
    ]

    molen=0
    SELFCOMP_MUMMER2BED ( input,molen )
}
