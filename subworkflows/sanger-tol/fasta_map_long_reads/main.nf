include { FASTXALIGN_PYFASTXINDEX    } from '../../../modules/sanger-tol/fastxalign/pyfastxindex/main'
include { FASTXALIGN_MINIMAP2ALIGN   } from '../../../modules/sanger-tol/fastxalign/minimap2align/main'
include { MINIMAP2_INDEX             } from '../../../modules/nf-core/minimap2/index/main'

include { BAM_SAMTOOLS_MERGE_MARKDUP } from '../bam_samtools_merge_markdup/main'

workflow FASTA_MAP_LONG_READS {

    take:
    ch_assemblies             // Channel [meta, assembly]
    ch_fasta                  // Channel [meta, fasta] OR [meta, [fasta1, fasta2, ..., fasta_n]]
    val_reads_per_fasta_chunk // integer: Number of reads per FASTA chunk for mapping
    val_output_bam            // boolean: if true output alignments in BAM format

    main:
    ch_versions = Channel.empty()

    //
    // Module: Index FASTA files
    //
    FASTXALIGN_PYFASTXINDEX(ch_fasta.transpose())
    ch_versions = ch_versions.mix(FASTXALIGN_PYFASTXINDEX.out.versions)

    //
    // Logic: Identify FASTA chunks
    //
    ch_fastx_chunks = FASTXALIGN_PYFASTXINDEX.out.index
        | map { meta, index, count ->
            def intcount = count.toInteger()
            def size     = val_reads_per_fasta_chunk
            def n_bins   = intcount.intdiv(size)
            chunkn       = (0..n_bins).collect()
            slices       = chunkn.collect { chunk ->
                def lower = chunk * size
                def upper = [lower + size - 1, intcount].min()

                return [ lower, upper ]
            }

            return [ meta, index, chunkn, slices ]
        }

    //
    // Logic: Count the total number of cram chunks for downstream grouping
    //
    ch_n_fasta_chunks = ch_fastx_chunks
        | map { meta, index, chunkn, _slices -> [ meta, chunkn ] }
        | transpose()
        | groupTuple(by: 0)
        | map { meta, chunkns -> [ meta, chunkns.size() ] }

    //
    // MODULE: generate minimap2 mmi file
    //
    MINIMAP2_INDEX(ch_assemblies)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    //
    // Module: Map slices of each FASTA file to the reference
    //
    ch_fasta_with_slices = ch_fasta
        | combine(ch_fastx_chunks, by: 0)
        | combine(MINIMAP2_INDEX.out.index, by: 0)
        | transpose()

    FASTXALIGN_MINIMAP2ALIGN(
        ch_fasta_with_slices,
        val_output_bam
    )
    ch_versions = ch_versions.mix(FASTXALIGN_MINIMAP2ALIGN.out.versions)

    ch_merge_input = FASTXALIGN_MINIMAP2ALIGN.out.bam
        | combine(ch_n_fasta_chunks, by: 0)
        | map { meta, bam, n_chunks ->
            def key = groupKey(meta, n_chunks)
            [key, bam]
        }
        | groupTuple(by: 0)
        | map { key, bam -> [key.target, bam] } // Get meta back out of groupKey

    //
    // Logic: Wrap this in the conditional so we don't unnecessarily run
    //        samtools faidx if no bam output
    //
    ch_output_bam = Channel.empty()

    if(val_output_bam) {
        BAM_SAMTOOLS_MERGE_MARKDUP(
            ch_merge_input,
            ch_assemblies,
            false
        )
        ch_versions = BAM_SAMTOOLS_MERGE_MARKDUP.out.versions

        ch_output_bam = ch_output_bam.mix(BAM_SAMTOOLS_MERGE_MARKDUP.out.bam)
    }

    emit:
    bam      = ch_output_bam
    paf      = FASTXALIGN_MINIMAP2ALIGN.out.paf
    versions = ch_versions
}
