#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURGEDUPS_FASTK } from '../../../../modules/purgedups/fastk/main.nf'

workflow test_purgedups_fastk {
    
    data = [
        [ id:'test'], // meta map
        [
           file('https://tolit.cog.sanger.ac.uk/test-data/Diarsia_rubi/working/ilDiaRubi1.hicanu.20220110/subset/ilDiaRubi1.trimmedReads.fasta.gz', checkIfExists: true),
        ],
    ]

    kmer=31
    myoutdir='/lustre/scratch123/tol/teams/grit/yy5/nf_test_small'
    FASTKDB='FASTKDB'

    PURGEDUPS_FASTK ( data,kmer,myoutdir,FASTKDB )
}
