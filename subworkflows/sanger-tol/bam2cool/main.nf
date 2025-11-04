include { SAMTOOLS_MARKDUP          } from '../../../modules/nf-core/samtools/markdup/main'
include { BAMTOBEDSORT            } from '../../../modules/sanger-tol/bam2bedsort/main'
include { CONTACTBED              } from '../../../modules/sanger-tol/contactbed/main'
include { GENERATE_CONTACTS_INDEX } from '../../../modules/sanger-tol/generatecontactsindex/main'
include { COOLER_CLOAD            } from '../../../modules/nf-core/cooler/cload/main'
include { COOLER_MERGE            } from '../../../modules/nf-core/cooler/merge/main'
include { COOLER_ZOOMIFY          } from '../../../modules/nf-core/cooler/zoomify/main'

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
    // Generate index file from contacts
    //
    GENERATE_CONTACTS_INDEX(
        CONTACTBED.out.bed
    )
    ch_versions = ch_versions.mix(GENERATE_CONTACTS_INDEX.out.versions)
    
    //
    // Generate individual .cool files
    //
    COOLER_CLOAD(
        GENERATE_CONTACTS_INDEX.out.contacts_with_index,
        ch_chrom_sizes,
        'full',
        val_bin_size
    )
    ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)

    //
    // Collect all individual .cool files for merging
    //
    ch_cool_files_for_merge = COOLER_CLOAD.out.cool
        .collect()
        .map { cool_files ->
            // Extract meta and paths from tuples
            def sample_ids = []
            def cool_paths = []
            cool_files.each { item ->
                if (item instanceof List && item.size() >= 2) {
                    if (item[0] != null && item[0].id != null) {
                        sample_ids.add(item[0].id)
                    }
                    if (item[1] != null) {
                        cool_paths.add(item[1])
                    }
                }
            }
            def merged_meta = [
                id: 'merged',
                samples: sample_ids
            ]
            [merged_meta, cool_paths]
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
