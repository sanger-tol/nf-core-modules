include { SAMTOOLS_MARKDUP          } from '../../../modules/nf-core/samtools/markdup/main'
include { BEDTOOLS_BAMTOBEDSORT     } from '../../../modules/nf-corel/bedtools/bamtobedsort/main'
include { GET_PAIRED_CONTACT_BED    } from '../../../modules/sanger-tol/hic/getpairedcontactbed/main'
include { COOLER_CLOAD              } from '../../../modules/sanger-tol/cooler/cload/main'
include { COOLER_MERGE              } from '../../../modules/sanger-tol/cooler/merge/main'
include { COOLER_ZOOMIFY            } from '../../../modules/sanger-tol/cooler/zoomify/main'

workflow BAM2COOL {

    take:
    ch_bam_list         // channel: [ val(meta), [ bam1, bam2, ... ] ]
    ch_reference        // channel: [ val(meta), path(fasta) ]
    ch_chrom_sizes      // channel: [ val(meta), path(chrom_sizes) ]
    val_bin_size        // value: bin size for cooler (e.g., 1000)

    main:

    ch_versions = Channel.empty()

    //
    // LOGIC: Transpose BAM list to process each BAM individually
    //
    ch_individual_bams = ch_bam_list
        .transpose()
        .map { meta, bam ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_${bam.baseName}"
            [new_meta, bam]
        }

    //
    // MODULE: SAMTOOLS MARKDUP - Mark duplicates in each BAM file
    //
    SAMTOOLS_MARKDUP(
        ch_individual_bams,
        ch_reference
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions)

    //
    // MODULE: BAMTOBED AND SORT - Convert each markdup BAM to BED individually
    //
    BEDTOOLS_BAMTOBEDSORT(
        SAMTOOLS_MARKDUP.out.bam
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBEDSORT.out.versions)

    //
    // MODULE: GENERATE PAIRED CONTACT BED - Extract Hi-C contact pairs
    //
    GET_PAIRED_CONTACT_BED(
        BEDTOOLS_BAMTOBEDSORT.out.sorted_bed
    )
    ch_versions = ch_versions.mix(GET_PAIRED_CONTACT_BED.out.versions)

    //
    // LOGIC: PREPARE COOLER INPUT - Combine contacts with chromosome sizes and bin size
    //
    ch_cooler_input = GET_PAIRED_CONTACT_BED.out.paired_contacts_bed
        .combine(ch_chrom_sizes)
        .map { meta_bam, contacts, meta_chrom, chrom_sizes ->
            [meta_bam, contacts, chrom_sizes, val_bin_size]
        }

    //
    // MODULE: COOLER CLOAD - Generate individual cool files
    //
    COOLER_CLOAD(
        ch_cooler_input
    )
    ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)

    //
    // LOGIC: COLLECT ALL COOL FILES FOR MERGING
    //
    ch_cool_files_for_merge = COOLER_CLOAD.out.cool
        .map { meta, cool -> cool }
        .collect()
        .map { cool_files ->
            def meta = [id: 'merged']
            [meta, cool_files]
        }

    //
    // MODULE: COOLER MERGE - Merge all individual cool files
    //
    COOLER_MERGE(
        ch_cool_files_for_merge
    )
    ch_versions = ch_versions.mix(COOLER_MERGE.out.versions)

    //
    // MODULE: COOLER ZOOMIFY - Create multi-resolution mcool file
    //
    COOLER_ZOOMIFY(
        COOLER_MERGE.out.cool
    )
    ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

    emit:
    individual_cool = COOLER_CLOAD.out.cool        // channel: [ val(meta), path(cool) ]
    merged_cool     = COOLER_MERGE.out.cool        // channel: [ val(meta), path(cool) ]
    mcool           = COOLER_ZOOMIFY.out.mcool     // channel: [ val(meta), path(mcool) ]
    versions        = ch_versions                  // channel: [ versions.yml ]
}
