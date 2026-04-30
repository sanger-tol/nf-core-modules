include { MINIMAP2_ALIGN     } from '../../../modules/sanger-tol/minimap2/align/main'
include { FIND_CONCATENATE   } from '../../../modules/nf-core/find/concatenate/main'
include { BEDTOOLS_SORT      } from '../../../modules/nf-core/bedtools/sort/main'
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

    // 2. Run patched minimap2 and emit one BED per read file
    MINIMAP2_ALIGN ( ch_reads_for_align, ch_reference, false, [], false, false, true )
    ch_paf_bed = MINIMAP2_ALIGN.out.bed

    // Group per-sample PAF/BED outputs into lists for downstream concatenation
    ch_paf_bed_grouped = ch_paf_bed
        .groupTuple(by: 0)
        .map { meta, paf_bed_files -> [ meta, paf_bed_files.sort { f -> f.name } ] }

    // 3. Merge all per-read BED outputs into one BED per sample
    FIND_CONCATENATE(ch_paf_bed_grouped)
    ch_versions = ch_versions.mix(FIND_CONCATENATE.out.versions_find)
    ch_versions = ch_versions.mix(FIND_CONCATENATE.out.versions_pigz)
    ch_versions = ch_versions.mix(FIND_CONCATENATE.out.versions_coreutils)

    BEDTOOLS_SORT(FIND_CONCATENATE.out.file_out, [])
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions_bedtools)

    // 4. Generate BedGraph from merged, sorted BED
    // Then call the module
    BEDTOOLS_GENOMECOV (
        BEDTOOLS_SORT.out.sorted.map { meta, bed -> [meta, bed, 1] },
        ch_chromsizes,
        'bedgraph',
        true
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions_bedtools)

    // 5. Convert to BigWig (Compulsory)
    UCSC_BEDGRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.genomecov,
        ch_chromsizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions_ucsc)

    emit:
    bigwig         = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
    bedgraph       = save_bedgraph ? BEDTOOLS_GENOMECOV.out.genomecov : Channel.empty()
    versions       = ch_versions
}
