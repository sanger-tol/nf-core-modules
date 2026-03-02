//
// Extract the masked regions from a Fasta file as BED,
// and prepare indexes for it
//

include { MASK_SOFTMASK2BED } from '../../../modules/sanger-tol/mask/softmask2bed/main'
include { BGZIPTABIX        } from '../../../modules/sanger-tol/bgziptabix/main'

workflow SOFT_MASKED_FASTA_REPEATS {
    take:
    ch_fasta_sequence_length // Channel [ val(meta), path(fasta), val(max_seq_length) ]

    main:

    ch_fasta = ch_fasta_sequence_length.map { meta, fasta, _max_seq_length -> [meta, fasta] }
    ch_bed = MASK_SOFTMASK2BED(ch_fasta).bed

    ch_bed_with_seq_length = ch_bed
        .join(ch_fasta_sequence_length)
        .map { meta, bed, _fasta, max_seq_length -> [meta, bed, max_seq_length] }
    BGZIPTABIX(ch_bed_with_seq_length)

    ch_repeats = BGZIPTABIX.out.gz_index
        .join(BGZIPTABIX.out.tbi, by: 0, remainder: true)
        .join(BGZIPTABIX.out.csi, by: 0, remainder: true)

    emit:
    repeats = ch_repeats // channel: [ meta, bed.gz, bed.gz.gzi, tbi?, csi? ]
}
