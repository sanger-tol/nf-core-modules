include { MINIMAP2_ALIGN     } from '../../../modules/nf-core/minimap2/align/main'
include { CAT_CAT            } from '../../../modules/nf-core/cat/cat/main'
include { BEDTOOLS_GENOMECOV } from '../../../modules/nf-core/bedtools/genomecov/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow READ_COVERAGE {

    take:
    ch_reads      // channel: [ val(meta), [ path(reads) ] ]
    ch_reference  // channel: [ val(meta2), path(fasta) ]
    ch_chromsizes // channel: [ path(chromsizes) ]
    save_bedgraph // boolean: true/false (Default: false)

    main:
    ch_versions = Channel.empty()

    // 1. Normalise read input to one tuple per read file
    //    Accept both:
    //      - [meta, path(read)]
    //      - [meta, [path(read1), path(read2), ...]]
    ch_reads_for_align = ch_reads
        .map { meta, reads ->
            def read_list = (reads instanceof List) ? reads : [reads]
            tuple(meta, read_list)
        }
        .flatMap { meta, read_list ->
            read_list.collect { read_file -> tuple(meta, read_file) }
        }

    // 2. Run Patched Minimap2 (outputs one PAF/BED per read file)
    MINIMAP2_ALIGN ( ch_reads_for_align, ch_reference, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_paf_bed = MINIMAP2_ALIGN.out.paf

    // Group per-sample PAF/BED outputs into lists for downstream consumers
    ch_paf_bed_grouped = ch_paf_bed
        .groupTuple(by: 0)
        .map { meta, paf_bed_files -> [ meta, paf_bed_files.sort { f -> f.name } ] }

    // 3. Merge all per-read BED outputs into one BED per sample
    CAT_CAT(ch_paf_bed_grouped)

    // 4. Generate BedGraph from merged BED
    BEDTOOLS_GENOMECOV ( 
        CAT_CAT.out.file_out,
        ch_chromsizes, 
        'bedgraph' 
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

    // 5. Convert to BigWig (Compulsory)
    UCSC_BEDGRAPHTOBIGWIG ( 
        BEDTOOLS_GENOMECOV.out.genomecov, 
        ch_chromsizes 
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions)

    emit:
    paf_bed        = ch_paf_bed                                          // Channel [meta, paf_or_bed] per input read file
    paf_bed_grouped = ch_paf_bed_grouped                                 // Channel [meta, [paf_or_bed1, ...]] per sample
    merged_bed     = CAT_CAT.out.file_out                                // Channel [meta, merged_bed]
    bigwig         = UCSC_BEDGRAPHTOBIGWIG.out.bw
    bedgraph       = save_bedgraph ? BEDTOOLS_GENOMECOV.out.genomecov : Channel.empty()
    versions       = ch_versions
}
