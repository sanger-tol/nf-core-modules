//
// MODULE IMPORT BLOCK
//
include { FINDTELOMERE   } from '../../../modules/sanger-tol/telomere/findtelomere/main'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'


workflow TELO_FINDER {

    take:
    ch_reference        // Channel [ val(meta), path(fasta) ]
    ch_telomereseq      // Channel [ val(meta), val(telomere_motif) ]
    val_split_telomere  // bool
    zip_bed             // bool — bgzip + tabix strand *.telomere.bed and windows outputs

    main:

    ch_joined = ch_reference
        .combine(ch_telomereseq, by: 0)
        .map { meta, reference, telomereseq -> tuple(meta, reference, telomereseq) }

    FINDTELOMERE(ch_joined, val_split_telomere)

    ch_windows_for_zip = val_split_telomere
        ? FINDTELOMERE.out.windows_fwd.mix(FINDTELOMERE.out.windows_rev)
        : FINDTELOMERE.out.windows

    if (zip_bed) {
        ch_beds_windows_for_zip_raw = FINDTELOMERE.out.telomere_bed_fwd
            .mix(FINDTELOMERE.out.telomere_bed_rev)
            .mix(ch_windows_for_zip)

        /*
         * Optional `path` outputs can arrive as an empty list when the glob matches nothing.
         * TABIX_BGZIPTABIX expects a single `path`, so normalise to (meta, path) tuples only.
         */
        ch_beds_windows_for_zip = ch_beds_windows_for_zip_raw
            .flatMap { meta, item ->
                def items = (item instanceof List) ? item : [ item ]
                items
                    .findAll { it != null }
                    .collect { file -> tuple(meta, file) }
            }

        TABIX_BGZIPTABIX(ch_beds_windows_for_zip)
    }

    emit:
    telomere         = FINDTELOMERE.out.telomere
    telomere_bed_fwd = FINDTELOMERE.out.telomere_bed_fwd
    telomere_bed_rev = FINDTELOMERE.out.telomere_bed_rev
    windows          = val_split_telomere ? Channel.empty()              : FINDTELOMERE.out.windows
    windows_fwd      = val_split_telomere ? FINDTELOMERE.out.windows_fwd : Channel.empty()
    windows_rev      = val_split_telomere ? FINDTELOMERE.out.windows_rev : Channel.empty()
    gz_index         = zip_bed ? TABIX_BGZIPTABIX.out.gz_index : Channel.empty()

}
