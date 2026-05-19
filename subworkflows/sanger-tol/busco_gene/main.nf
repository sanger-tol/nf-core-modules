//
// Run BUSCO in genome mode on a FASTA assembly, then split gene coordinates from
// full_table.tsv into Complete / Duplicated / Fragmented bedGraph-style tables (sequence, start, end, score).
//
include { BUSCO_BUSCO                      } from '../../../modules/nf-core/busco/busco/main'
include { BUSCOFULLTABLETOGENEBEDGRAPH     } from '../../../modules/sanger-tol/busco/buscofulltabletogenebedgraph/main'
include { TABIX_BGZIPTABIX                 } from '../../../modules/nf-core/tabix/bgziptabix/main'


workflow BUSCO_GENE {

    take:
    ch_fasta                  // channel: [ val(meta), path(fasta) ]
    ch_busco_lineage          // channel: [ val(meta), val(lineage_string) ]
    val_busco_lineages_path   // value: path to local BUSCO lineage DB (--download_path). Use [] to download.
    val_zip_bedgraph          // bool — bgzip + tabix the three *.bedgraph outputs

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

    BUSCOFULLTABLETOGENEBEDGRAPH(BUSCO_BUSCO.out.full_table)

    if (val_zip_bedgraph) {
        ch_bedgraphs_for_zip_raw = BUSCOFULLTABLETOGENEBEDGRAPH.out.complete_bedgraph
            .mix(BUSCOFULLTABLETOGENEBEDGRAPH.out.duplicated_bedgraph)
            .mix(BUSCOFULLTABLETOGENEBEDGRAPH.out.fragmented_bedgraph)

        ch_bedgraphs_for_zip = ch_bedgraphs_for_zip_raw
            .flatMap { meta, item ->
                def items = (item instanceof List) ? item : [ item ]
                items.collect { file -> tuple(meta, file) }
            }

        TABIX_BGZIPTABIX(ch_bedgraphs_for_zip)
    }

    emit:
    complete_bedgraph      = BUSCOFULLTABLETOGENEBEDGRAPH.out.complete_bedgraph
    duplicated_bedgraph    = BUSCOFULLTABLETOGENEBEDGRAPH.out.duplicated_bedgraph
    fragmented_bedgraph    = BUSCOFULLTABLETOGENEBEDGRAPH.out.fragmented_bedgraph
    full_table             = BUSCO_BUSCO.out.full_table
    busco_dir              = BUSCO_BUSCO.out.busco_dir
    batch_summary          = BUSCO_BUSCO.out.batch_summary
    short_summaries_txt    = BUSCO_BUSCO.out.short_summaries_txt
    short_summaries_json   = BUSCO_BUSCO.out.short_summaries_json
    busco_log              = BUSCO_BUSCO.out.log
    gz_index               = val_zip_bedgraph ? TABIX_BGZIPTABIX.out.gz_index : Channel.empty()
}
