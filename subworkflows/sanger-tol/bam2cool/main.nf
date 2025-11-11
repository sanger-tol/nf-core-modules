include { BEDTOOLS_BAMTOBEDSORT     } from '../../../modules/sanger-tol/bedtools/bamtobedsort/main'
include { CONTACTBED                } from '../../../modules/sanger-tol/contactbed/main'
include { GENERATE_CONTACTS_INDEX   } from '../../../modules/sanger-tol/generatecontactsindex/main'
include { COOLER_CLOAD              } from '../../../modules/nf-core/cooler/cload/main'
include { COOLER_MERGE              } from '../../../modules/nf-core/cooler/merge/main'
include { COOLER_ZOOMIFY            } from '../../../modules/nf-core/cooler/zoomify/main'

workflow BAM2COOL {

    take:
    ch_bam_list     // channel: [val(meta), [bam1, bam2, ...]]
    ch_chrom_sizes // channel: [val(meta), path(chrom_sizes)]
    val_bin_size    // integer: bin size for cooler

    main:

    ch_versions = Channel.empty()

    //
    // Convert marked BAMs to sorted BED
    //
    BEDTOOLS_BAMTOBEDSORT(
        ch_bam_list
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBEDSORT.out.versions)

    //
    // Generate paired contacts BED
    //
    CONTACTBED(
        BEDTOOLS_BAMTOBEDSORT.out.sorted_bed
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
        .groupTuple(by: 0)
        .map { cool_files ->
            def cool_paths = []
            cool_files.each { item ->
                if (item instanceof List && item.size() >= 2 && item[1] != null) {
                    cool_paths.add(item[1])
                }
            }
            def merged_meta = [ id: 'merged' ]
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
