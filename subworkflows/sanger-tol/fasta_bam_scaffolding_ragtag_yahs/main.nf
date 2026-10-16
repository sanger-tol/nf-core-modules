include { BEDTOOLS_BAMTOBEDSORT                      } from '../../../modules/sanger-tol/bedtools/bamtobedsort/main'
include { RAGTAG_SCAFFOLD                            } from '../../../modules/nf-core/ragtag/scaffold/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_CONTIGS   } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SCAFFOLDS } from '../../../modules/nf-core/samtools/faidx/main'
include { YAHS as YAHS_PRIMARY_HAP                   } from '../../../modules/nf-core/yahs/main'
include { YAHS as YAHS_RAGTAG_HAPS                   } from '../../../modules/nf-core/yahs/main'
include { YAHS_MAKEPAIRSFILE                         } from '../../../modules/sanger-tol/yahs/makepairsfile/main'

workflow FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS {
    take:
    ch_fasta // [meta,  [hap1, hap2, ..., hap_n]]
    ch_hic_bam // [meta, [hap1_bam, hap2_bam, ..., hap_n_bam]]

    main:
    ch_input = ch_fasta
        .combine(ch_hic_bam, by: 0)
        .flatMap { meta, fastas, bams ->
            [fastas, bams]
                .transpose()
                .withIndex()
                .collect { hap_bam, index ->
                    def out_meta = meta + [_hap: "hap${index + 1}"]
                    [out_meta, hap_bam[0], hap_bam[1]]
                }
        }
        .multiMap { meta, fasta, bam ->
            assemblies: [meta, fasta]
            bam: [meta, bam]
        }

    //
    // Module: Convert BAM to name-sorted BED
    //
    BEDTOOLS_BAMTOBEDSORT(ch_input.bam)

    //
    // Module: Index input assemblies
    //
    SAMTOOLS_FAIDX_CONTIGS(
        ch_input.assemblies.map { meta, assembly -> [meta, assembly, []] },
        false,
    )

    //
    // Module: scaffold primary assembly contigs with YaHS
    //
    ch_yahs_primary_input = ch_input.assemblies
        .filter { meta, _asm -> meta._hap == "hap1" }
        .combine(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        .combine(BEDTOOLS_BAMTOBEDSORT.out.sorted_bed, by: 0)
        .map { meta, fasta, fai, bed ->
            // add empty AGP input
            [meta, fasta, fai, bed, []]
        }

    YAHS_PRIMARY_HAP(ch_yahs_primary_input)

    //
    // Module: scaffold additional haplotypes against primary haplotype
    //
    ch_ragtag_inputs = YAHS_PRIMARY_HAP.out.scaffolds_fasta
        .combine(ch_input.assemblies.filter { meta, _asm -> meta._hap != "hap1" })
        .filter { meta_primary, _primary_scaffolds, meta_hap, _hap_contigs ->
            (meta_primary - meta_primary.subMap("_hap")) == meta_hap - meta_hap.subMap("_hap")
        }
        .multiMap { meta_primary, primary_scaffolds, meta_hap, hap_contigs ->
            ragtag_reference: [meta_primary, primary_scaffolds]
            ragtag_contigs: [meta_hap, hap_contigs]
        }

    RAGTAG_SCAFFOLD(
        ch_ragtag_inputs.ragtag_contigs,
        ch_ragtag_inputs.ragtag_reference,
        [[], []],
        [[], [], []],
    )

    //
    // Module: run YaHS on the RagTag scaffolds
    //
    ch_yahs_ragtag_inputs = ch_input.assemblies
        .filter { meta, _asm -> meta._hap != "hap1" }
        .combine(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        .combine(BEDTOOLS_BAMTOBEDSORT.out.sorted_bed, by: 0)
        .combine(RAGTAG_SCAFFOLD.out.corrected_agp, by: 0)

    YAHS_RAGTAG_HAPS(ch_yahs_ragtag_inputs)

    //
    // Logic: Re-combine all YaHS outputs into a single set of channels
    //
    ch_yahs_scaffolds = YAHS_PRIMARY_HAP.out.scaffolds_fasta.mix(YAHS_RAGTAG_HAPS.out.scaffolds_fasta)
    ch_yahs_agp = YAHS_PRIMARY_HAP.out.scaffolds_agp.mix(YAHS_RAGTAG_HAPS.out.scaffolds_agp)
    ch_yahs_initial_break_agp = YAHS_PRIMARY_HAP.out.initial_break_agp.mix(YAHS_RAGTAG_HAPS.out.initial_break_agp)
    ch_yahs_round_agp = YAHS_PRIMARY_HAP.out.round_agp.mix(YAHS_RAGTAG_HAPS.out.round_agp)
    ch_yahs_binary = YAHS_PRIMARY_HAP.out.binary.mix(YAHS_RAGTAG_HAPS.out.binary)
    ch_yahs_log = YAHS_PRIMARY_HAP.out.log.mix(YAHS_RAGTAG_HAPS.out.log)

    //
    // Module: Index output scaffolds
    //
    SAMTOOLS_FAIDX_SCAFFOLDS(
        ch_yahs_scaffolds.map { meta, assembly -> [meta, assembly, []] },
        true,
    )

    //
    // Module: Make pairs file to build contact maps with
    //
    ch_pairs_input = SAMTOOLS_FAIDX_SCAFFOLDS.out.fai
        .combine(ch_yahs_agp, by: 0)
        .combine(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        .combine(BEDTOOLS_BAMTOBEDSORT.out.sorted_bed, by: 0)

    YAHS_MAKEPAIRSFILE(ch_pairs_input)

    emit:
    scaffolds_fasta      = ch_yahs_scaffolds
    scaffolds_agp        = ch_yahs_agp
    scaffolds_pairs      = YAHS_MAKEPAIRSFILE.out.pairs
    scaffolds_chromsizes = SAMTOOLS_FAIDX_SCAFFOLDS.out.sizes
    ragtag_scaffolds     = RAGTAG_SCAFFOLD.out.corrected_assembly
    ragtag_agp           = RAGTAG_SCAFFOLD.out.corrected_agp
    ragtag_stats         = RAGTAG_SCAFFOLD.out.corrected_stats
    ragtag_paf           = RAGTAG_SCAFFOLD.out.ragtag_paf
    ragtag_delta         = RAGTAG_SCAFFOLD.out.ragtag_delta
    confidence_txt       = RAGTAG_SCAFFOLD.out.confidence_txt
    yahs_bin             = ch_yahs_binary
    yahs_inital          = ch_yahs_initial_break_agp
    yahs_intermediate    = ch_yahs_round_agp
    yahs_log             = ch_yahs_log
}
