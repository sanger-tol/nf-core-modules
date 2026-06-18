//
// Run BUSCO in genome mode on a FASTA assembly, then split gene coordinates from
// full_table.tsv into Complete / Duplicated / Fragmented bedGraph-style tables (sequence, start, end, score).
//
include { BUSCO_BUSCO                           } from '../../../modules/nf-core/busco/busco/main'
include { BUSCO_BUSCOFULLTABLETOGENEBEDGRAPH    } from '../../../modules/sanger-tol/busco/buscofulltabletogenebedgraph/main'
include { HTSLIB_BGZIPTABIX                     } from '../../../modules/nf-core/htslib/bgziptabix/main'


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

    BUSCO_BUSCOFULLTABLETOGENEBEDGRAPH(BUSCO_BUSCO.out.full_table)

    if (val_zip_bedgraph) {
        ch_bedgraphs_for_zip_raw = BUSCO_BUSCOFULLTABLETOGENEBEDGRAPH.out.complete_bedgraph
            .mix(BUSCO_BUSCOFULLTABLETOGENEBEDGRAPH.out.duplicated_bedgraph)
            .mix(BUSCO_BUSCOFULLTABLETOGENEBEDGRAPH.out.fragmented_bedgraph)

        ch_bedgraphs_for_zip = ch_bedgraphs_for_zip_raw
            .flatMap { meta, item ->
                def items = (item instanceof List) ? item : [ item ]
                items.collect { file -> tuple(meta, file, [], []) }
            }

        HTSLIB_BGZIPTABIX(
            ch_bedgraphs_for_zip,
            'compress',
            true,
            'bedgraph'
        )

        ch_gz_index = HTSLIB_BGZIPTABIX.out.output
            .combine(HTSLIB_BGZIPTABIX.out.index)
            .filter { meta, gz, _meta2, idx -> idx.name.startsWith(gz.name) }
            .map { meta, gz, _meta2, idx -> tuple(meta, gz, idx) }

    } else {
        ch_gz_index = channel.empty()
    }

    emit:
    complete_bedgraph      = BUSCO_BUSCOFULLTABLETOGENEBEDGRAPH.out.complete_bedgraph
    duplicated_bedgraph    = BUSCO_BUSCOFULLTABLETOGENEBEDGRAPH.out.duplicated_bedgraph
    fragmented_bedgraph    = BUSCO_BUSCOFULLTABLETOGENEBEDGRAPH.out.fragmented_bedgraph
    full_table             = BUSCO_BUSCO.out.full_table
    busco_dir              = BUSCO_BUSCO.out.busco_dir
    batch_summary          = BUSCO_BUSCO.out.batch_summary
    short_summaries_txt    = BUSCO_BUSCO.out.short_summaries_txt
    short_summaries_json   = BUSCO_BUSCO.out.short_summaries_json
    busco_log              = BUSCO_BUSCO.out.log
    gz_index               = ch_gz_index
}
