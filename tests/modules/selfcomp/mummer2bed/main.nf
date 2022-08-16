#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SELFCOMP_MUMMER2BED } from '../../../../modules/selfcomp/mummer2bed/main.nf'

workflow test_selfcomp_mummer2bed {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://tolit.cog.sanger.ac.uk/test-data/Oscheius_sp/genomic_data/nxOscSpen2/selfcomp/nxOscSpen2.mummer", checkIfExists: true)
    ]

    molen=0
    SELFCOMP_MUMMER2BED ( input,molen )
}
