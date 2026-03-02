//
// Prepare all the indexes for a Fasta file
//

include { SAMTOOLS_BGZIP } from '../../../modules/nf-core/samtools/bgzip/main'
include { SAMTOOLS_DICT  } from '../../../modules/nf-core/samtools/dict/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'

workflow FASTA_COMPRESS_INDEX {
    take:
    ch_fasta // Channel [ val(meta), path(fasta) ]
    val_run_samtools_dict // Boolean

    main:

    // Compress the Fasta file
    ch_compressed_fasta = SAMTOOLS_BGZIP(ch_fasta).fasta

    // Generate .gzi and .fai index files, and chromosome size file
    ch_fasta_with_dummy_fai = ch_compressed_fasta.map { meta, fasta -> [meta, fasta, []] }
    ch_samtools_faidx = SAMTOOLS_FAIDX(ch_fasta_with_dummy_fai, true).fai

    // Read the .fai file, extract sequence statistics, and make an extended meta map
    ch_sequence_map = ch_samtools_faidx.map { meta, fai ->
        [meta, meta + get_sequence_map(fai)]
    }
    // Update all channels to use the extended meta map
    ch_fasta_gz = ch_compressed_fasta.join(ch_sequence_map).map { _meta, path, ext_meta -> [ext_meta, path] }
    ch_faidx = ch_samtools_faidx.join(ch_sequence_map).map { _meta, path, ext_meta -> [ext_meta, path] }
    ch_gzi = SAMTOOLS_FAIDX.out.gzi.join(ch_sequence_map).map { _meta, path, ext_meta -> [ext_meta, path] }
    ch_sizes = SAMTOOLS_FAIDX.out.sizes.join(ch_sequence_map).map { _meta, path, ext_meta -> [ext_meta, path] }

    // Generate Samtools dictionary
    ch_samtools_dict = channel.empty()
    if (val_run_samtools_dict) {
        ch_samtools_dict = SAMTOOLS_DICT(ch_fasta_gz).dict
    }

    emit:
    fasta_gz = ch_fasta_gz // Channel [ val(meta), path(genome.fa.gz) ]
    faidx    = ch_faidx // Channel [ val(meta), path(genome.fa.gz.fai) ]
    dict     = ch_samtools_dict // Channel [ val(meta), path(genome.fa.gz.dict) ]
    gzi      = ch_gzi // Channel [ val(meta), path(genome.fa.gz.gzi) ]
    sizes    = ch_sizes // Channel [ val(meta), path(genome.fa.gz.sizes) ]
}

// Read the .fai file to extract the number of sequences, the maximum and total sequence length
// Inspired from https://github.com/nf-core/rnaseq/blob/3.10.1/lib/WorkflowRnaseq.groovy
def get_sequence_map(fai_file) {
    def n_sequences = 0
    def max_length = 0
    def total_length = 0
    fai_file.eachLine { line ->
        def lspl = line.split('\t')
        // def chrom  = lspl[0]
        def length = lspl[1].toLong()
        n_sequences += 1
        total_length += length
        if (length > max_length) {
            max_length = length
        }
    }

    def sequence_map = [:]
    sequence_map.n_sequences = n_sequences
    sequence_map.total_length = total_length
    if (n_sequences) {
        sequence_map.max_length = max_length
    }
    return sequence_map
}
