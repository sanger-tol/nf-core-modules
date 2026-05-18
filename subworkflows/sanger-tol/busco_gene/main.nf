//
// Run BUSCO in genome mode on a FASTA assembly, then split gene coordinates from
// full_table.tsv into Complete / Duplicated / Fragmented bedGraph-style tables (sequence, start, end, score).
//
include { BUSCO_BUSCO                         } from '../../../modules/nf-core/busco/busco/main'
include { BUSCO_FULLTABLE_TO_GENE_BEDGRAPH } from '../../../modules/sanger-tol/busco/fulltabletogenebedgraph/main'


workflow BUSCO_GENE {

    take:
    ch_fasta                  // channel: [ val(meta), path(fasta) ]
    ch_busco_lineage          // channel: [ val(meta), val(lineage_string) ]
    val_busco_lineages_path   // value: path to local BUSCO lineage DB (--download_path). Use [] to download.

    main:

    ch_busco_inputs = ch_fasta
        .join(ch_busco_lineage, by: 0)
        .multiMap { meta, fasta, lineage ->
            assemblies: tuple(meta, fasta)
            lineage: lineage
        }

    BUSCO_BUSCO(
        ch_busco_inputs.assemblies,
        'genome',
        ch_busco_inputs.lineage,
        val_busco_lineages_path ?: [],
        [],
        true
    )

    BUSCO_FULLTABLE_TO_GENE_BEDGRAPH(BUSCO_BUSCO.out.full_table)

    emit:
    complete_bedgraph      = BUSCO_FULLTABLE_TO_GENE_BEDGRAPH.out.complete_bedgraph      // channel: [ val(meta), path(complete_buscos.bedgraph) ]
    duplicated_bedgraph    = BUSCO_FULLTABLE_TO_GENE_BEDGRAPH.out.duplicated_bedgraph    // channel: [ val(meta), path(duplicated_buscos.bedgraph) ]
    fragmented_bedgraph    = BUSCO_FULLTABLE_TO_GENE_BEDGRAPH.out.fragmented_bedgraph    // channel: [ val(meta), path(fragmented_buscos.bedgraph) ]
    full_table             = BUSCO_BUSCO.out.full_table                                   // channel: [ val(meta), path(full_table.tsv) ]
    busco_dir              = BUSCO_BUSCO.out.busco_dir
    batch_summary          = BUSCO_BUSCO.out.batch_summary
    short_summaries_txt    = BUSCO_BUSCO.out.short_summaries_txt
    short_summaries_json   = BUSCO_BUSCO.out.short_summaries_json
    busco_log              = BUSCO_BUSCO.out.log
}
