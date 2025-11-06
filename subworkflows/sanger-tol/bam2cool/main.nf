include { SAMTOOLS_MARKDUP } from '../../../modules/nf-core/samtools/markdup/main'
include { BAMTOBEDSORT     } from '../../../modules/sanger-tol/bam2bedsort/main'
include { CONTACTBED       } from '../../../modules/sanger-tol/contactbed/main'
include { COOLER_CLOAD     } from '../../../modules/nf-core/cooler/cload/main'
include { COOLER_MERGE     } from '../../../modules/nf-core/cooler/merge/main'
include { COOLER_ZOOMIFY   } from '../../../modules/nf-core/cooler/zoomify/main'

workflow BAM2COOL {

    take:
    ch_bam_list     // channel: [val(meta), [bam1, bam2, ...]]
    ch_reference    // channel: [val(meta), path(fasta)]
    ch_chrom_sizes // channel: [val(meta), path(chrom_sizes)]
    val_bin_size    // integer: bin size for cooler

    main:

    ch_versions = Channel.empty()

    //
    //  Mark duplicates for each BAM
    //
    SAMTOOLS_MARKDUP(
        ch_bam_list,
        ch_reference
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions)

    SAMTOOLS_MARKDUP.out.bam.view()
    //
    // Convert marked BAMs to sorted BED
    //
    BAMTOBEDSORT(
        SAMTOOLS_MARKDUP.out.bam
    )
    ch_versions = ch_versions.mix(BAMTOBEDSORT.out.versions)

    //
    // Generate paired contacts BED
    //
    CONTACTBED(
        BAMTOBEDSORT.out.sorted_bed
    )
    ch_versions = ch_versions.mix(CONTACTBED.out.versions)

    //
    // Prepare input for individual COOLER_CLOAD
    //
    ch_cooler_input = CONTACTBED.out.paired_contacts_bed.combine(ch_chrom_sizes)
        .map { meta_bam, contacts, meta_chrom, chrom_sizes_path ->
            [meta_bam, contacts, chrom_sizes_path, val_bin_size]
        }

    //
    // Generate individual .cool files
    //
    COOLER_CLOAD(
        ch_cooler_input
    )
    ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)

    //
    // Collect all individual .cool files for merging
    //
    ch_cool_files_for_merge = COOLER_CLOAD.out.cool
        .collect()
        .map { cool_files ->
            // Preserve list of contributing sample IDs
            def merged_meta = [
                id: 'merged',
                samples: cool_files*.first.id
            ]
            // Extract only the cool file paths
            [merged_meta, cool_files*.last]
        }

    //
    // Merge individual .cool files
    //
    COOLER_MERGE(
        ch_cool_files_for_merge
    )
    ch_versions = ch_versions.mix(COOLER_MERGE.out.versions)

    //
    // Generate multi-resolution .mcool
    //
    COOLER_ZOOMIFY(
        COOLER_MERGE.out.cool
    )
    ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

    //
    // Emit final outputs
    //
    emit:
    merged_cool     = COOLER_MERGE.out.cool
    mcool           = COOLER_ZOOMIFY.out.mcool
    versions        = ch_versions
}
