#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { MINIPROT_INDEX } from '../../../../modules/miniprot/index/main.nf'
include { MINIPROT_ALIGN } from '../../../../modules/miniprot/align/main.nf'

workflow test_miniprot_align {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) 
    input_pep = file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true)

    MINIPROT_INDEX ( [ [id:'test'], input] )

    input_ref = MINIPROT_INDEX.out.index

    paf_format = false
    gff_format = true
    gtf_format = false
    
    MINIPROT_ALIGN ( input_ref, input_pep, paf_format, gff_format, gtf_format )
}
