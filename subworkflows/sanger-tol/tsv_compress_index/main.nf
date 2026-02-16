include { TABIX_BGZIP                    } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as TABIX_TABIX_CSI } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TABIX_TBI } from '../../../modules/nf-core/tabix/tabix/main'

workflow TSV_COMPRESS_INDEX {
    take:
    ch_tsv // channel: [ val(meta), [ tsv ] ]
    ch_max_seq_length // channel: [ val(meta), val(max_seq_length) ]

    main:

    // Compress the TSV file
    ch_compressed_tsv = TABIX_BGZIP(ch_tsv).output

    // Try indexing the TSV file in two formats for maximum compatibility
    // but each has its own limitations
    tabix_selector = ch_compressed_tsv
        .join(ch_max_seq_length, by: 0, remainder: true)
        .branch { meta, tsv, max_seq_length ->
            lonely_max_length: tsv == null
            [meta, max_seq_length]
            no_index: max_seq_length == null || max_seq_length >= 2 ** 32
            [meta, tsv]
            tbi_and_csi: max_seq_length < 2 ** 29
            [meta, tsv]
            only_csi: true
            // 2**29 <= max_seq_length < 2**32
            [meta, tsv]
        }

    // Output channels to tell the downstream subworkflows which indexes are missing
    no_csi = tabix_selector.no_index
    no_tbi = tabix_selector.no_index.mix(tabix_selector.only_csi)

    // Do the indexing on the compatible TSV files
    for_csi = tabix_selector.tbi_and_csi.mix(tabix_selector.only_csi)
    for_tbi = tabix_selector.tbi_and_csi
    ch_csi = TABIX_TABIX_CSI(for_csi).index
    ch_tbi = TABIX_TABIX_TBI(for_tbi).index

    // Combine the compressed TSV with the indexes for downstream subworkflows
    ch_combined = ch_compressed_tsv
        .join(ch_tbi, by: 0, remainder: true)
        .join(ch_csi, by: 0, remainder: true)

    emit:
    gz       = ch_compressed_tsv // channel: [ val(meta), [ tsv.gz ] ]
    tbi      = ch_tbi // channel: [ val(meta), [ tsv.gz.tbi ] ]
    csi      = ch_csi // channel: [ val(meta), [ tsv.gz.csi ] ]
    no_tbi   = no_tbi // channel: [ val(meta), [ tsv.gz ] ]
    no_csi   = no_csi // channel: [ val(meta), [ tsv.gz ] ]
    combined = ch_combined // channel: [ tsv.gz, tbi?, csi? ]
}
