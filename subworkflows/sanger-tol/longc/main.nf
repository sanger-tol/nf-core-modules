//
// Subworkflow: Long-C alignment (optional digest + minimap2 + annotate + pairs + contact maps)
//

include { LONGC_DIGESTREADS                } from '../../../modules/sanger-tol/longc/digestreads/main'
include { LONGC_ANNOTATEFRAG               } from '../../../modules/sanger-tol/longc/annotatefrag/main'
include { LONGC_PAIRTOOLSPARSE2           } from '../../../modules/sanger-tol/longc/pairtoolsparse2/main'
include { FIND_CONCATENATE                 } from '../../../modules/sanger-tol/find/concatenate/main'
include { PAIRTOOLS_RESTRICT               } from '../../../modules/nf-core/pairtools/restrict/main'
include { PAIRTOOLS_SORT                   } from '../../../modules/nf-core/pairtools/sort/main'
include { PAIRTOOLS_DEDUP                  } from '../../../modules/nf-core/pairtools/dedup/main'
include { MINIMAP2_INDEX                   } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                   } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_FAIDX                   } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MERGE                   } from '../../../modules/nf-core/samtools/merge/main'
include { COOLER_DIGEST                    } from '../../../modules/nf-core/cooler/digest/main'
include { COOLER_CLOAD                     } from '../../../modules/nf-core/cooler/cload/main'
include { COOLER_ZOOMIFY                   } from '../../../modules/nf-core/cooler/zoomify/main'
include { PRETEXTMAP                       } from '../../../modules/nf-core/pretextmap/main'

workflow LONGC {

    take:
    reference           // channel: [ val(meta), path(fasta) ]
    longc_reads         // channel: [ val(meta), path(reads) ] or [ val(meta), [ path(reads), ... ] ]
    val_skip_digest     // bool — align raw reads with original read names
    val_cutter          // string, e.g. 'NlaIII'
    val_restrict_frags  // bool — restrict pairs to enzyme fragments (requires digest)
    val_cool            // bool — build .cool contact map
    val_mcool           // bool — build multi-resolution .mcool (requires val_cool)
    val_pretext         // bool — build Pretext contact map
    val_cool_bin_size   // int — cooler bin size in bp
    val_dedup           // bool — sort and deduplicate pairs (recommended before contact maps)

    main:

    //
    // Chromosome sizes for pairtools parse2 / cooler (also used when merging CRAM/BAM lists)
    //
    SAMTOOLS_FAIDX(
        reference.map { meta, fasta -> [ meta, fasta, [] ] },
        true
    )

    ch_reads_list = longc_reads.map { meta, reads ->
        def reads_list = (reads instanceof List) ? reads.sort { a, b -> a.name <=> b.name } : [ reads ]
        tuple(meta, reads_list)
    }

    //
    // Optionally digest reads (split concatemers at restriction sites)
    // skip_digest=true: align raw reads directly (merge multi-file lists first)
    // skip_digest=false: digest_reads.py on FASTQ, or BAM/CRAM via samtools fastq
    //
    if (val_skip_digest) {
        ch_reads_list.branch { meta, reads_list ->
            single: reads_list.size() == 1
                return tuple(meta, reads_list[0])
            multi_fastq: reads_list.every { f -> f.name.toLowerCase() ==~ /.*\.(fq|fastq|fa|fasta)(\.gz)?$/ }
                return tuple(meta, reads_list)
            multi_bam: reads_list.every { f -> f.name.toLowerCase().endsWith('.bam') }
                return tuple(meta, reads_list)
            multi_cram: reads_list.every { f -> f.name.toLowerCase().endsWith('.cram') }
                return tuple(meta, reads_list)
            invalid: true
                error("LONGC: reads list for sample '${meta.id}' mixes file types or uses unsupported extensions")
        }
        .set { ch_reads_branch }

        FIND_CONCATENATE(
            ch_reads_branch.multi_fastq,
            false
        )

        ch_fai_gzi = SAMTOOLS_FAIDX.out.fai
            .join(SAMTOOLS_FAIDX.out.gzi, by: 0, remainder: true)
            .map { meta, fai, gzi -> tuple(meta, fai, gzi ?: []) }

        ch_bam_merge_input = ch_reads_branch.multi_bam
            .mix(ch_reads_branch.multi_cram)
            .combine(reference, by: 0)
            .combine(ch_fai_gzi, by: 0)
            .multiMap { meta, bams, fasta, fai, gzi ->
                bam:   tuple(meta, bams, [])
                fasta: tuple(meta, fasta, fai, gzi)
            }

        SAMTOOLS_MERGE(
            ch_bam_merge_input.bam,
            ch_bam_merge_input.fasta
        )

        ch_reads_for_align = ch_reads_branch.single
            .mix(FIND_CONCATENATE.out.file_out)
            .mix(SAMTOOLS_MERGE.out.bam)
    }
    else {
        ch_for_digest = ch_reads_list.join(reference, by: 0)

        LONGC_DIGESTREADS(ch_for_digest, val_cutter)
        ch_reads_for_align = LONGC_DIGESTREADS.out.digested_reads
        ch_digest_versions = LONGC_DIGESTREADS.out.versions_python
            .mix(LONGC_DIGESTREADS.out.versions_digest_reads)
            .mix(LONGC_DIGESTREADS.out.versions_samtools)
    }

    //
    // Index reference FASTA with minimap2, then align reads to the index
    //
    MINIMAP2_INDEX(reference)

    ch_for_align = ch_reads_for_align
        .combine(MINIMAP2_INDEX.out.index, by: 0)

    MINIMAP2_ALIGN(
        ch_for_align.map { meta, reads, index -> [ meta, reads ] },
        ch_for_align.map { meta, reads, index -> [ meta, index ] },
        true,    // bam_format
        'bai',   // bam_index_extension
        false,   // cigar_paf_format
        false    // cigar_bam
    )

    //
    // Annotate aligned fragments: group by read, filter, order, assign fragment IDs
    //
    ch_align_bam = MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index, by: 0)

    LONGC_ANNOTATEFRAG(ch_align_bam)

    //
    // Convert BAM to pairs (pairtools parse2)
    //
    ch_bam_for_parse = LONGC_ANNOTATEFRAG.out.bam

    LONGC_PAIRTOOLSPARSE2(
        ch_bam_for_parse.combine(SAMTOOLS_FAIDX.out.sizes, by: 0)
    )

    //
    // Optionally restrict pairs to restriction fragments (cooler digest → pairtools restrict)
    //
    use_restrict = !val_skip_digest && val_restrict_frags

    if (use_restrict) {
        ch_digest = reference.combine(SAMTOOLS_FAIDX.out.sizes, by: 0)

        COOLER_DIGEST(
            ch_digest.map { meta, fasta, sizes -> fasta },
            ch_digest.map { meta, fasta, sizes -> sizes },
            val_cutter
        )

        ch_restriction_bed = ch_digest
            .map { meta, fasta, sizes ->
                def bed_name = "${fasta.baseName}_${val_cutter.replaceAll(/[^0-9a-zA-Z]+/, '_')}.bed"
                tuple(bed_name, meta)
            }
            .join(COOLER_DIGEST.out.bed.map { bed -> tuple(bed.name, bed) })
            .map { bed_name, meta, bed -> tuple(meta, bed) }

        ch_pairs_bed = LONGC_PAIRTOOLSPARSE2.out.pairs
            .join(ch_restriction_bed, by: 0)

        PAIRTOOLS_RESTRICT(
            ch_pairs_bed.map { meta, pairs, bed -> [ meta, pairs ] },
            ch_pairs_bed.map { meta, pairs, bed -> bed }
        )

        ch_pairs_final = PAIRTOOLS_RESTRICT.out.restrict
    }
    else {
        ch_pairs_final = LONGC_PAIRTOOLSPARSE2.out.pairs
    }

    //
    // Optionally deduplicate pairs (pairtools sort → dedup; wf-pore-c hi_c path)
    //
    if (val_dedup) {
        PAIRTOOLS_SORT(ch_pairs_final)
        PAIRTOOLS_DEDUP(PAIRTOOLS_SORT.out.sorted)
        ch_pairs_out = PAIRTOOLS_DEDUP.out.pairs
    }
    else {
        ch_pairs_out = ch_pairs_final
    }

    //
    // Convert pairs to cool / mcool / pretext (optional)
    //
    use_cool = val_cool || val_mcool

    if (use_cool) {
        ch_cooler_input = ch_pairs_out
            .combine(SAMTOOLS_FAIDX.out.sizes, by: 0)
            .multiMap { meta, pairs, sizes ->
                pairs: [ meta, pairs, [] ]
                sizes: [ meta, sizes ]
            }

        COOLER_CLOAD(
            ch_cooler_input.pairs,
            ch_cooler_input.sizes,
            'pairs',
            val_cool_bin_size
        )
    }

    if (val_mcool && val_cool) {
        COOLER_ZOOMIFY(COOLER_CLOAD.out.cool)
    }

    if (val_pretext) {
        PRETEXTMAP(
            ch_pairs_out,
            [[], [], []]
        )
    }

    ch_versions = Channel.empty()
        .mix(MINIMAP2_INDEX.out.versions_minimap2)
        .mix(MINIMAP2_ALIGN.out.versions_minimap2)
        .mix(LONGC_ANNOTATEFRAG.out.versions_annotate_frag)
        .mix(LONGC_PAIRTOOLSPARSE2.out.versions_pairtools)
    if (val_skip_digest) {
        ch_versions = ch_versions
            .mix(FIND_CONCATENATE.out.versions_find)
            .mix(FIND_CONCATENATE.out.versions_pigz)
            .mix(FIND_CONCATENATE.out.versions_coreutils)
            .mix(SAMTOOLS_MERGE.out.versions_samtools)
    }
    if (!val_skip_digest) {
        ch_versions = ch_versions.mix(ch_digest_versions)
    }
    if (use_restrict) {
        ch_versions = ch_versions
            .mix(COOLER_DIGEST.out.versions_cooler)
            .mix(PAIRTOOLS_RESTRICT.out.versions_pairtools)
    }
    if (val_dedup) {
        ch_versions = ch_versions
            .mix(PAIRTOOLS_SORT.out.versions_pairtools)
            .mix(PAIRTOOLS_DEDUP.out.versions_pairtools)
    }
    if (use_cool) {
        ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions_cooler)
    }
    if (val_mcool && val_cool) {
        ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions_cooler)
    }
    if (val_pretext) {
        ch_versions = ch_versions.mix(PRETEXTMAP.out.versions_pretextmap)
    }

    emit:
    bam      = LONGC_ANNOTATEFRAG.out.bam.join(LONGC_ANNOTATEFRAG.out.index, by: 0)
    pairs    = ch_pairs_out
    cool     = use_cool ? COOLER_CLOAD.out.cool : Channel.empty()
    mcool    = (val_mcool && val_cool) ? COOLER_ZOOMIFY.out.mcool : Channel.empty()
    pretext  = val_pretext ? PRETEXTMAP.out.pretext : Channel.empty()
    versions = ch_versions
}
