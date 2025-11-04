include { BEDTOOLS_BAMTOBEDSORT                      } from '../../../modules/sanger-tol/bedtools/bamtobedsort/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_CONTIGS   } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SCAFFOLDS } from '../../../modules/nf-core/samtools/faidx/main'
include { YAHS                                       } from '../../../modules/nf-core/yahs/main'
include { YAHS_MAKEPAIRSFILE                         } from '../../../modules/sanger-tol/yahs/makepairsfile/main'

include { PAIRS_CREATE_CONTACT_MAPS                  } from '../pairs_create_contact_maps/main'

workflow FASTA_BAM_SCAFFOLDING_YAHS {
    take:
    ch_fasta      // [meta, assembly]
    ch_hic_bam    // [meta, bam/bed]
    val_cool_bin  // val: cooler cload parameter

    main:
    ch_versions = Channel.empty()

    //
    // Module: Convert BAM to name-sorted BED
    //
    BEDTOOLS_BAMTOBEDSORT(ch_hic_bam)
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBEDSORT.out.versions)

    //
    // Module: Index input assemblies
    //
    SAMTOOLS_FAIDX_CONTIGS(
        ch_fasta,
        [[],[]],
        false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_CONTIGS.out.versions)

    //
    // Module: scaffold contigs with YaHS
    //
    ch_yahs_input = ch_fasta
        | combine(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        | combine(ch_bed, by: 0)

    YAHS(ch_yahs_input)
    ch_versions = ch_versions.mix(YAHS.out.versions)

    //
    // Module: Index output scaffolds
    //
    SAMTOOLS_FAIDX_SCAFFOLDS(
        YAHS.out.scaffolds_fasta,
        [[],[]],
        true
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SCAFFOLDS.out.versions)

    //
    // Module: Make pairs file to build contact maps with
    //
    ch_pairs_input = SAMTOOLS_FAIDX_SCAFFOLDS.out.fai
        | join(YAHS.out.scaffolds_agp, by: 0)
        | join(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        | join(YAHS.out.binary, by: 0)

    YAHS_MAKEPAIRSFILE(ch_pairs_input)
    ch_versions = ch_versions.mix(YAHS_MAKEPAIRSFILE.out.versions)

    //
    // Subworkflow: Create Hi-C contact maps for visualisation of scaffolding outputs
    //
    PAIRS_CREATE_CONTACT_MAPS(
        YAHS_MAKEPAIRSFILE.out.pairs,
        SAMTOOLS_FAIDX_SCAFFOLDS.out.sizes,
        val_cool_bin
    )
    ch_versions = ch_versions.mix(PAIRS_CREATE_CONTACT_MAPS.out.versions)

    emit:
    assemblies  = YAHS.out.scaffolds_fasta
    agp         = YAHS.out.scaffolds_agp
    pretext     = PRETEXTMAP.out.pretext
    pretext_png = PRETEXTSNAPSHOT.out.image
    cool        = COOLER_ZOOMIFY.out.mcool
    hic         = JUICERTOOLS_PRE.out.hic
    versions    = ch_versions
}
