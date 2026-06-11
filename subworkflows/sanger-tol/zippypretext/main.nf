//
// MODULE IMPORT BLOCK
//
include { PRETEXT2ASM       } from '../../../modules/sanger-tol/pretext/pretext2asm/main'
include { JUICERC           } from '../../../modules/sanger-tol/pretext/juicerc/main'
include { PRETEXT_MAKEPAIRS } from '../../../modules/sanger-tol/pretext/makepairs/main'
include { PRETEXTMAP        } from '../../../modules/nf-core/pretextmap/main'

workflow ZIPPYPRETEXT {

    take:
    ch_fasta   // channel: [ val(meta), path(fasta) ]
    ch_agp     // channel: [ val(meta), path(agp) ] — PretextView-curated AGP
    ch_hicmap  // channel: [ val(meta), path(hicmap) ] — Hi-C `.bin` or BAM
    ch_idxfile // channel: [ val(meta), path(fai) ] — contig-level FASTA index

    main:

    ch_zippy_inputs = ch_fasta
        .combine(ch_agp, by: 0)
        .combine(ch_hicmap, by: 0)
        .combine(ch_idxfile, by: 0)
        .multiMap { meta, fasta, agp, hicmap, idx ->
            fasta:   tuple(meta, fasta)
            agp:     tuple(meta, agp)
            hicmap:  tuple(meta, hicmap)
            idxfile: tuple(meta, idx)
        }

    PRETEXT_PRETEXT2ASM (
        ch_zippy_inputs.fasta,
        ch_zippy_inputs.agp.map { meta, agp -> agp }
    )

    PRETEXT_JUICERC (
        ch_zippy_inputs.hicmap,
        PRETEXT_PRETEXT2ASM.out.correctedagp.map { meta, agp -> agp },
        ch_zippy_inputs.idxfile.map { meta, idx -> idx }
    )

    ch_makepairs_input = PRETEXT_JUICERC.out.alignment
        .combine(PRETEXT_JUICERC.out.outlog)
        .map { meta, alignment, outlog -> tuple(meta, alignment, outlog) }

    PRETEXT_MAKEPAIRS(ch_makepairs_input)

    PRETEXTMAP(
        PRETEXT_MAKEPAIRS.out.pairs,
        [[], [], []]
    )

    emit:
    pairs        = PRETEXT_MAKEPAIRS.out.pairs
    pretext      = PRETEXTMAP.out.pretext

}
