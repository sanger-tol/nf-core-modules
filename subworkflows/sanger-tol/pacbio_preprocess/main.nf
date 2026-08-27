include { BLAST_BLASTN                         } from '../../../modules/sanger-tol/blast/blastn/main'
include { BLAST_MAKEBLASTDB                    } from '../../../modules/nf-core/blast/makeblastdb/main'
include { HIFITRIMMER_PROCESSBLAST             } from '../../../modules/nf-core/hifitrimmer/processblast/main'
include { HIFITRIMMER_TRIM                     } from '../../../modules/nf-core/hifitrimmer/trim/main'
include { LIMA as LIMA_ULI                     } from '../../../modules/nf-core/lima/main'
include { LIMA as LIMA_PIMMS                   } from '../../../modules/nf-core/lima/main'
include { PBMARKDUP                            } from '../../../modules/nf-core/pbmarkdup/main'
include { UNTAR                                } from '../../../modules/nf-core/untar/main'

workflow PACBIO_PREPROCESS {

    take:
    ch_reads_standard     // Channel [meta, input]: standard reads → optional hifi-trimmer
    ch_reads_uli          // Channel [meta, input]: ULI reads → LIMA(global) + pbmarkdup
    ch_reads_pimms        // Channel [meta, input]: pimms reads → per-sample LIMA + pbmarkdup
    ch_reads_amplified    // Channel [meta, input]: amplified reads → pbmarkdup only
    ch_adapter_yaml       // Channel [meta, yaml]: yaml for hifitrimmer (applies to ch_reads_standard)
    ch_pimms_adapters     // Channel [meta, path]: per-sample LIMA adapter file for pimms
    val_hifi_adapter      // Path to Hifi adapter DB or Hifi adapter fasta to make database for blastn
    val_uli_adapter       // Adapter file for ULI LIMA

    main:
    lima_reports    = channel.empty()
    lima_summary    = channel.empty()
    pbmarkdup_stats = channel.empty()

    //
    // ULI: LIMA with global adapter
    //
    LIMA_ULI(ch_reads_uli, val_uli_adapter)
    lima_reports = lima_reports.mix(LIMA_ULI.out.report)
    lima_summary = lima_summary.mix(LIMA_ULI.out.summary)
    ch_uli_post_lima = LIMA_ULI.out.bam
        .mix(LIMA_ULI.out.fastq)
        .mix(LIMA_ULI.out.fasta)
        .mix(LIMA_ULI.out.fastqgz)
        .mix(LIMA_ULI.out.fastagz)

    //
    // PiMmS: per-sample LIMA with matching adapter file
    //
    ch_pimms_lima_input = ch_reads_pimms
        .join(ch_pimms_adapters, by: 0, remainder: true)
        .map { meta, reads, adapter ->
            if (!adapter) { error "PiMmS adapter file is required for ${meta.id}: ${reads}" }
            [ meta, reads, adapter ]
        }
        .multiMap { meta, reads, adapter ->
            reads:    [ meta, reads ]
            adapters: adapter
        }

    LIMA_PIMMS(ch_pimms_lima_input.reads, ch_pimms_lima_input.adapters)
    lima_reports = lima_reports.mix(LIMA_PIMMS.out.report)
    lima_summary = lima_summary.mix(LIMA_PIMMS.out.summary)
    ch_pimms_post_lima = LIMA_PIMMS.out.bam
        .mix(LIMA_PIMMS.out.fastq)
        .mix(LIMA_PIMMS.out.fasta)
        .mix(LIMA_PIMMS.out.fastqgz)
        .mix(LIMA_PIMMS.out.fastagz)

    //
    // pbmarkdup: ULI (post-LIMA) + pimms (post-LIMA) + amplified (raw)
    //
    PBMARKDUP(ch_uli_post_lima.mix(ch_pimms_post_lima).mix(ch_reads_amplified))
    pbmarkdup_stats = pbmarkdup_stats.mix(PBMARKDUP.out.log)

    //
    // TRIMMING WITH HIFITRIMMER
    //
    hifitrimmer_summary = channel.empty()
    hifitrimmer_bed     = channel.empty()
    trimmed_bam         = channel.empty()
    trimmed_cram        = channel.empty()
    trimmed_sam         = channel.empty()
    trimmed_fasta       = channel.empty()
    trimmed_fastq       = channel.empty()
    ch_hifitrimmer_input = ch_reads_standard.mix(PBMARKDUP.out.markduped)

    if ( val_hifi_adapter ) {

        ch_input_skip_trim = ch_hifitrimmer_input
            .join(ch_adapter_yaml, by: 0, remainder: true)
            .filter { _meta, _reads, yaml -> !yaml }
            .map { meta, reads, _yaml -> [meta, reads] }

        ch_input_skip_trim
            .subscribe { _meta, reads ->
                log.warn "No adapter YAML provided, skipping adapter trimming step for: ${reads}"
            }

        ch_input_to_trim = ch_hifitrimmer_input
            .combine(ch_adapter_yaml, by: 0)
            .map { meta, reads, _yaml -> [meta, reads] }

        adapter_fasta_ch = channel.of([ [id: file(val_hifi_adapter).baseName], file(val_hifi_adapter) ])
        if ( val_hifi_adapter.endsWith('.tar.gz') ) {
            UNTAR( adapter_fasta_ch )
            adapter_db = UNTAR.out.untar
        } else {
            BLAST_MAKEBLASTDB( adapter_fasta_ch, [] )
            adapter_db = BLAST_MAKEBLASTDB.out.db
        }

        BLAST_BLASTN ( ch_input_to_trim, adapter_db.collect(), [],[],[] )

        ch_input_processblast = BLAST_BLASTN.out.txtgz.combine( ch_adapter_yaml, by: 0 )
            .multiMap { meta, blastn, yaml ->
                blastn: [ meta, blastn ]
                yaml: [ meta, yaml ]
            }

        HIFITRIMMER_PROCESSBLAST ( ch_input_processblast.blastn, ch_input_processblast.yaml )

        hifitrimmer_summary = hifitrimmer_summary.mix ( HIFITRIMMER_PROCESSBLAST.out.summary )
        hifitrimmer_bed     = hifitrimmer_bed.mix ( HIFITRIMMER_PROCESSBLAST.out.bed )

        ch_input_filterbam = ch_input_to_trim.combine( HIFITRIMMER_PROCESSBLAST.out.bed, by: 0 )
        HIFITRIMMER_TRIM ( ch_input_filterbam )

        trimmed_bam   = trimmed_bam.mix( HIFITRIMMER_TRIM.out.bam )
        trimmed_cram  = trimmed_cram.mix( HIFITRIMMER_TRIM.out.cram )
        trimmed_sam   = trimmed_sam.mix( HIFITRIMMER_TRIM.out.sam )
        trimmed_fasta = trimmed_fasta.mix( HIFITRIMMER_TRIM.out.fasta )
        trimmed_fastq = trimmed_fastq.mix( HIFITRIMMER_TRIM.out.fastq )
    } else {
        ch_input_skip_trim = ch_hifitrimmer_input
    }

    ch_input_skip_trim_branch = ch_input_skip_trim
        .branch { meta, reads ->
            bam: reads.name.endsWith('.bam')
                return [ meta, reads ]
            fastx: true
                return [ meta, reads ]
        }

    emit:
    untrimmed_fastx     = ch_input_skip_trim_branch.fastx   // [meta, fastx] untrimmed reads in FASTA/FASTQ format
    untrimmed_bam       = ch_input_skip_trim_branch.bam     // [meta, bam] untrimmed reads in BAM format
    trimmed_cram        = trimmed_cram                      // [meta, CRAM] preprocessed reads in CRAM format
    trimmed_bam         = trimmed_bam                       // [meta, BAM] preprocessed reads in BAM format
    trimmed_sam         = trimmed_sam                       // [meta, SAM] preprocessed reads in SAM format
    trimmed_fasta       = trimmed_fasta                     // [meta, FASTA] preprocessed reads in FASTA format
    trimmed_fastq       = trimmed_fastq                     // [meta, FASTQ] preprocessed reads in FASTQ format
    lima_report         = lima_reports                      // [meta, report]
    lima_summary        = lima_summary                      // [meta, summary]
    hifitrimmer_bed     = hifitrimmer_bed                   // [meta, bed]
    hifitrimmer_summary = hifitrimmer_summary               // [meta, summary]
    pbmarkdup_stats     = pbmarkdup_stats                   // [meta, log]
}
