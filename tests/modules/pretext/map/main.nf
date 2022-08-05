#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRETEXT_MAP } from '../../../../modules/pretext/map/main.nf'

workflow test_pretext_map {
    
/*    input = [
        [[ id:'test' ], // meta map
        file(params.tol_test_data['test']['pEimTen1']['assembly']['alignments_sorted'], checkIfExists: true)],
        file(params.tol_test_data['test']['pEimTen1']['assembly']['scaffolds_fai'], checkIfExists: true)
    ]
*/

input = [
    [ id:'test' ],
    file("/lustre/scratch124/tol/projects/darwin/data/protists/Eimeria_tenella/working/pEimTen1.scaff/out.break.salsa/alignments_sorted.txt")]

fai = file("/lustre/scratch124/tol/projects/darwin/data/protists/Eimeria_tenella/working/pEimTen1.scaff/out.nobreak.salsa/scaffolds_FINAL.fasta.fai")
    PRETEXT_MAP ( input , fai)
}
