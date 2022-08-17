#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { SELFCOMP_MAPIDS } from '../../../../modules/selfcomp/mapids/main.nf'
workflow test_selfcomp_mapids {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.tol_test_data['small_genome']['Oscheius_sp']['genomic_data']['nxOscSpen2']['selfcomp']['mummer_bed'], checkIfExists: true)
    ]
    agp = file(params.tol_test_data['small_genome']['Oscheius_sp']['genomic_data']['nxOscSpen2']['selfcomp']['split_agp'], checkIfExists: true)
    SELFCOMP_MAPIDS ( input, agp )
}
