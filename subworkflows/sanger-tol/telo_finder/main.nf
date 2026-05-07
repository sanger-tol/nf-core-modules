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

    FINDTELOMERE(
        ch_joined.map{ meta, fasta, motif -> [meta, fasta] }, 
        ch_joined.map{ meta, fasta, motif -> motif },         
        val_split_telomere                                    
        )

    ch_windows_for_zip = val_split_telomere
        ? FINDTELOMERE.out.windows_fwd.mix(FINDTELOMERE.out.windows_rev)
        : FINDTELOMERE.out.windows

    ch_beds_windows_for_zip = FINDTELOMERE.out.telomere_bed_fwd
        .mix(FINDTELOMERE.out.telomere_bed_rev)
        .mix(ch_windows_for_zip)

    TABIX_BGZIPTABIX(
        ch_beds_windows_for_zip.filter { _meta, _file -> zip_bed }
    )

    emit:
    telomere         = FINDTELOMERE.out.telomere
    telomere_bed_fwd = FINDTELOMERE.out.telomere_bed_fwd
    telomere_bed_rev = FINDTELOMERE.out.telomere_bed_rev
    windows          = val_split_telomere ? Channel.empty()              : FINDTELOMERE.out.windows
    windows_fwd      = val_split_telomere ? FINDTELOMERE.out.windows_fwd : Channel.empty()
    windows_rev      = val_split_telomere ? FINDTELOMERE.out.windows_rev : Channel.empty()
    gz_index         = TABIX_BGZIPTABIX.out.gz_index

}
