#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAKECMAP_RENAMECMAPIDS } from '../../../../modules/makecmap/renamecmapids/main.nf'

workflow test_makecmap_renamecmapids {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file("/lustre/scratch123/tol/resources/nextflow/test-data/Mucor_piriformis/genomic_data/gzMucPiri1/bionano/gzMucPiri1_DLE1.cmap", checkIfExists: true)
    ]

    keyfile = file("/lustre/scratch123/tol/resources/nextflow/test-data/Mucor_piriformis/genomic_data/gzMucPiri1/bionano/gzMucPiri1_DLE1_key.txt",checkIfExists: true)

    MAKECMAP_RENAMECMAPIDS ( input, keyfile)
}
