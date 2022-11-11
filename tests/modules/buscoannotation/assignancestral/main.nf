#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCOANNOTATION_ASSIGNANCESTRAL } from '../../../../modules/buscoannotation/assignancestral/main.nf'

workflow test_buscoannotation_assignancestral {
    
    ancestral_locations = [
        [ id:'test', single_end:false ], // meta map
        file("https://tolit.cog.sanger.ac.uk/test-data/Leptidea_sinapis/genomic_data/ilLepSina1/busco/buscopainter_complete_location.tsv",checkIfExists: true)
    ]

    fulltable = file("https://tolit.cog.sanger.ac.uk/test-data/Leptidea_sinapis/genomic_data/ilLepSina1/busco/ilLepSina1_full_table.tsv",checkIfExists: true)

    BUSCOANNOTATION_ASSIGNANCESTRAL ( ancestral_locations, fulltable  )
}
