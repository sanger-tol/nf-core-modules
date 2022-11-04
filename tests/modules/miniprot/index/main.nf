#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIPROT_INDEX } from '../../../../modules/miniprot/index/main.nf'

workflow test_miniprot_index {

    /*fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)*/
    fasta =file('https://tolit.cog.sanger.ac.uk/test-data/Oscheius_sp/assembly/nxOscSpef2.fasta', checkIfExists: true)

    MINIPROT_INDEX ( [ [id:'test'], fasta ] )
}
