include { CRAMALIGN_GENCRAMCHUNKS         } from '../../../modules/sanger-tol/cramalign/gencramchunks'
include { CRAMALIGN_MINIMAP2ALIGN         } from '../../../modules/sanger-tol/cramalign/minimap2align/main'
include { MINIMAP2_INDEX                  } from '../../../modules/nf-core/minimap2/index/main'
include { SAMTOOLS_FAIDX                  } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX                  } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                  } from '../../../modules/nf-core/samtools/merge/main'

workflow LONGREAD_MAPPING {

    take:
    ch_assemblies           // Channel [meta, assembly]
    ch_crams                // Channel [meta, cram] OR [meta, [cram1, cram2, ..., cram_n]]
    val_cram_chunk_size     // integer: Number of CRAM slices per chunk for mapping

    main:
    ch_versions = Channel.empty()

    //
    // Logic: check if CRAM files are accompanied by an index
    //        Get indexes, and index those that aren't
    //
    ch_cram_raw = ch_crams
        | transpose()
        | branch { meta, cram ->
            def cram_file = file(cram, checkIfExists: true)
            def index = cram + ".crai"
            have_index: file(index).exists()
                return [ meta, cram_file, file(index, checkIfExists: true) ]
            no_index: true
                return [ meta, cram_file ]
        }

    //
    // Module: Index CRAM files without indexes
    //
    SAMTOOLS_INDEX(ch_cram_raw.no_index)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_cram_indexed = ch_cram_raw.have_index
        | mix(
            ch_cram_raw.no_index.join(SAMTOOLS_INDEX.out.crai)
        )

    //
    // Module: Process the cram index files to determine how many
    //         chunks to split into for mapping
    //
    ch_cram_indexes = ch_cram_indexed
        | map { meta, cram, index -> [ meta, index ] }

    CRAMALIGN_GENCRAMCHUNKS(
        ch_cram_indexes,
        val_cram_chunk_size
    )
    ch_versions = ch_versions.mix(CRAMALIGN_GENCRAMCHUNKS.out.versions)

    //
    // Logic: Count the total number of cram chunks for downstream grouping
    //
    ch_n_cram_chunks = CRAMALIGN_GENCRAMCHUNKS.out.cram_slices
        | map { meta, chunkn, _slices -> [ meta, chunkn ] }
        | transpose()
        | groupTuple(by: 0)
        | map { meta, chunkns -> [ meta, chunkns.size() ] }

    //
    // Logic: Re-join the cram files and indexes to their chunk information
    //
    ch_cram_with_slices = ch_cram_indexed
        | combine(CRAMALIGN_GENCRAMCHUNKS.out.cram_slices, by: 0)

    //
    // MODULE: generate minimap2 mmi file
    //
    MINIMAP2_INDEX(ch_assemblies)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    ch_cram_chunks = ch_cram_with_slices
        | transpose()
        | combine(MINIMAP2_INDEX.out.index, by: 0)

    CRAMALIGN_MINIMAP2ALIGN(ch_cram_chunks)
    ch_versions = ch_versions.mix(CRAMALIGN_MINIMAP2ALIGN.out.versions)

    //
    // Module: Index assembly fastas
    //
    SAMTOOLS_FAIDX(
        ch_assemblies, // reference
        [ [:],[] ],    // fai
        false          // get sizes
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    //
    // Logic: create a channel with both fai and gzi for each assembly
    //        We do it here so we don't cause downstream issues with the
    //        remainder join
    //
    ch_fai_gzi = SAMTOOLS_FAIDX.out.fai
        | join(SAMTOOLS_FAIDX.out.gzi, by: 0, remainder: true)
        | map { meta, fai, gzi -> [ meta, fai, gzi ?: [] ] }


    //
    // Logic: Prepare input for merging bams.
    //        We use the ch_n_cram_chunks to set a groupKey so that
    //        we emit groups downstream ASAP once all bams have been made
    //
    ch_samtools_merge_input = CRAMALIGN_MINIMAP2ALIGN.out.bam
        | combine(ch_n_cram_chunks, by: 0)
        | map { meta, bam, n_chunks ->
            def key = groupKey(meta, n_chunks)
            [key, bam]
        }
        | groupTuple(by: 0)
        | map { key, bam -> [key.target, bam] } // Get meta back out of groupKey
        | combine(ch_assemblies, by: 0)
        | combine(ch_fai_gzi, by: 0)
        | multiMap { meta, bams, assembly, fai, gzi ->
            bam:   [ meta, bams ]
            fasta: [ meta, assembly ]
            fai:   [ meta, fai ]
            gzi:   [ meta, gzi ]
        }

    //
    // Module: Merge position-sorted bam files
    //
    SAMTOOLS_MERGE(
        ch_samtools_merge_input.bam,
        ch_samtools_merge_input.fasta,
        ch_samtools_merge_input.fai,
        ch_samtools_merge_input.gzi,
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    emit:
    bam      = SAMTOOLS_MERGE.out.bam
    versions = ch_versions
}
