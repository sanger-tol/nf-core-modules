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

    main:
    ch_versions = channel.empty()

    // 1. Normalise read input to one tuple per read file
    //    Accept both:
    //      - [meta, path(read)]
    //      - [meta, [path(read1), path(read2), ...]]
    // 1. Normalise reads and pair each read with its matching reference (by meta id)
    ch_reads
        .map { meta, reads ->
            def read_list = (reads instanceof List) ? reads : [reads]
            tuple(meta, read_list)
        }
        .flatMap { meta, read_list ->
            read_list.withIndex().collect { read_file, idx ->
                def meta_with_read_idx = meta + [read_idx: idx]
                tuple(meta_with_read_idx, read_file)
            }
        }
        .combine(ch_reference, by: 0)
        .multiMap { _sample_id, meta, read_file, meta2, reference ->
            reads: tuple(meta, read_file)
            reference: tuple(meta2, reference)
        }
        .set { ch_align_split }

    // 3. Run minimap2 and emit one BED per read file
    MINIMAP2_ALIGN ( ch_align_split.reads, ch_align_split.reference, false, [], false, false, true )

    // Group per-sample PAF/BED outputs into lists for downstream concatenation
    ch_paf_bed_grouped = MINIMAP2_ALIGN.out.bed
        .groupTuple(by: 0)
        .map { meta, paf_bed_files -> [ meta, paf_bed_files.sort { f -> f.name } ] }

    // 4. Merge all per-read BED outputs into one BED per sample
    FIND_CONCATENATE(ch_paf_bed_grouped)

    // 5. Sort bed
    // 5. Sort bed
    BEDTOOLS_SORT(FIND_CONCATENATE.out.file_out, [])

    // 5. Generate BedGraph from merged, sorted BED
    // Then call the module
    BEDTOOLS_GENOMECOV (
        BEDTOOLS_SORT.out.sorted.map { meta, bed -> [meta, bed, 1] },
        ch_chromsizes,
        'bedgraph',
        true
    )

    // 6. Convert to BigWig (Compulsory)
    UCSC_BEDGRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.genomecov,
        ch_chromsizes
    )

    emit:
    bigwig         = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
    bedgraph       = BEDTOOLS_GENOMECOV.out.genomecov
}
