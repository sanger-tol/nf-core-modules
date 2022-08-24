#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLAST_FILTERBLAST } from '../../../../modules/blast/filterblast/main.nf'

workflow test_blast_filterblast {
    
    input = [
        [id: 'test', type: 'rna'], // meta map
        file(params.tol_test_data['small_genome']['Gae_host']['genomic_data']['blast']['blast_tsv'], checkIfExists: true)
    ]

    BLAST_FILTERBLAST ( input )
}
