#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAKECMAP_CMAP2BED } from '../../../../modules/makecmap/cmap2bed/main.nf'

workflow test_makecmap_cmap2bed {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.tol_test_data['small_genome']['Mucor_piriformis']['genomic_data']['bionano_cmap'], checkIfExists: true)
    ]
    enzyme = "DLE1"
    MAKECMAP_CMAP2BED ( input, enzyme )
}
