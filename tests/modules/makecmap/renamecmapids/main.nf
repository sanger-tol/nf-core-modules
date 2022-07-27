#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAKECMAP_RENAMECMAPIDS } from '../../../../modules/makecmap/renamecmapids/main.nf'

workflow test_makecmap_renamecmapids {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.tol_test_data['small_genome']['Mucor_piriformis']['genomic_data']['bionano_cmap'], checkIfExists: true)
    ]

    keyfile = file(params.tol_test_data['small_genome']['Mucor_piriformis']['genomic_data']['bionano_key'], checkIfExists: true)

    MAKECMAP_RENAMECMAPIDS ( input, keyfile)
}
