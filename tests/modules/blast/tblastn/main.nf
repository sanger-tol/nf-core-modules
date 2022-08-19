#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLAST_TBLASTN } from '../../../../modules/blast/tblastn/main.nf'
include { BLAST_MAKEBLASTDB } from '../../../../modules/blast/makeblastdb/main.nf'

workflow test_blast_tblastn {

    input = [ file(params.tol_test_data['small_genome']['Oscheius_sp']['genomic_data']['nxOscSpen2']['selfcomp']['split_fasta'], checkIfExists: true) ]
    input_pep = [ file(params.tol_test_data['small_genome']['Gae_host']['genomic_data']['pep']['pep_seq'], checkIfExists: true) ]
    BLAST_MAKEBLASTDB ( input )
    BLAST_TBLASTN ( [ [id:'test'], input_pep ], BLAST_MAKEBLASTDB.out.db )
}
