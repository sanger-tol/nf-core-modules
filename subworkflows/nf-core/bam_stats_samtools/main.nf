//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../../modules/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/samtools/flagstat/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    bam_bai // channel: [ val(meta), [ bam/cram ], [bai/csi] ]
    fasta   // channel: [ fasta ]

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_STATS ( bam_bai, fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    SAMTOOLS_FLAGSTAT ( bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    SAMTOOLS_IDXSTATS ( bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    emit:
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
