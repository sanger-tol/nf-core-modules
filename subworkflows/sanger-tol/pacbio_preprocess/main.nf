include { BLAST_BLASTN                         } from '../../../modules/nf-core/blast/blastn/main'
include { BLAST_MAKEBLASTDB                    } from '../../../modules/nf-core/blast/makeblastdb/main'
include { HIFITRIMMER_PROCESSBLAST             } from '../../../modules/nf-core/hifitrimmer/processblast/main'
include { HIFITRIMMER_FILTERBAM                } from '../../../modules/nf-core/hifitrimmer/filterbam/main'
include { LIMA                                 } from '../../../modules/nf-core/lima/main'
include { PBMARKDUP                            } from '../../../modules/nf-core/pbmarkdup/main'
include { SAMTOOLS_FASTA                       } from '../../../modules/nf-core/samtools/fasta/main'
include { SAMTOOLS_FASTQ                       } from '../../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_IMPORT as FQ2BAM            } from '../../../modules/nf-core/samtools/import/main'
include { SAMTOOLS_IMPORT as FQ2CRAM_UNTRIM    } from '../../../modules/nf-core/samtools/import/main'
include { SAMTOOLS_IMPORT as FQ2CRAM_TRIM      } from '../../../modules/nf-core/samtools/import/main'
include { SAMTOOLS_VIEW                        } from '../../../modules/nf-core/samtools/view/main'
include { SEQKIT_FQ2FA                         } from '../../../modules/nf-core/seqkit/fq2fa/main'

workflow PACBIO_PREPROCESS {

    take:
    ch_reads                    // Channel [meta, input]: input reads in FASTA/FASTQ/BAM format, if trimming, only FASTQ/BAM
    ch_adapter_yaml             // Channel [meta, yaml]: yaml file for hifitrimmer adapter trimming
    val_adapter_fasta           // Adapter fasta to make database for blastn
    val_uli_primers             // Primer file for lima
    val_pbmarkdup               // Options to run pbmarkdup
    val_output_format           // Output format: 'cram' or 'fastx' ( default: 'fastx' )

    main:
    ch_versions = channel.empty()
    lima_reports = channel.empty()
    lima_summary = channel.empty()
    pbmarkdup_stats = channel.empty()
    untrimmed_cram = channel.empty()
    untrimmed_fastx = channel.empty()
    trimmed_fastx = channel.empty()
    trimmed_cram = channel.empty()

    //
    // DEMULTIPLEX WITH LIMA
    //
    if ( val_uli_primers ) {
        LIMA( ch_reads, val_uli_primers )
        ch_versions = ch_versions.mix( LIMA.out.versions )

        lima_reports = lima_reports.mix( LIMA.out.report )
        lima_summary = lima_summary.mix( LIMA.out.summary )

        // prepare input for markdup or trimming
        ch_input_for_md = LIMA.out.bam
            .mix(LIMA.out.fastq)
            .mix(LIMA.out.fasta)
            .mix(LIMA.out.fastqgz)
            .mix( LIMA.out.fastagz )
    } else {
        ch_input_for_md = ch_reads
    }

    //
    // MARKDUP WITH PBMARKDUP
    //
    if ( val_pbmarkdup ) {
        PBMARKDUP( ch_input_for_md )
        ch_versions = ch_versions.mix( PBMARKDUP.out.versions )
        pbmarkdup_stats = pbmarkdup_stats.mix( PBMARKDUP.out.log )

        ch_input_to_trim = PBMARKDUP.out.markduped
    } else {
        // If not running markdup, pass the input to trimming step
        ch_input_to_trim = ch_input_for_md
    }

    //
    // TRIMMING WITH HIFITRIMMER
    //
    hifitrimmer_summary = Channel.empty()
    hifitrimmer_bed = Channel.empty()
    if ( val_adapter_fasta ) {
        // Assign ch_input_skip_trimm to those without adapter yaml for trimming
        ch_input_skip_trim = ch_input_to_trim
        .join(ch_adapter_yaml, by: 0, remainder: true)
        .filter { meta, reads, yaml -> yaml == null }
        .map { meta, reads, yaml -> [meta, reads ] }

        // Warning for skip trimming
        ch_input_skip_trim
            .subscribe { meta, reads ->
                log.warn "No adapter YAML provided, skipping adapter trimming step for: ${reads}"
            }

        // PREPARE INPUT FOR TRIMMING
        // Combine adapter yaml to input reads, only those with adapter yaml will be used for trimming, skip those without
        ch_input_trim_yaml = ch_input_to_trim
            .combine(ch_adapter_yaml, by: 0)
            .map { meta, reads, yaml -> [meta, reads] }

        // Prepare reads for filter
        ch_input_trim_branch = ch_input_trim_yaml
            .branch { meta, reads -> def filename = reads.toString()
                fasta: filename =~ /\.(fasta|fa|fna)(\.gz)?$/
                    return [ meta, reads ]
                fastq: filename =~ /\.(fastq|fq)(\.gz)?$/
                    return [ meta, reads ]
                bam: filename.endsWith('.bam')
                    return [ meta, reads ]
            }
        // Make adapter database
        BLAST_MAKEBLASTDB( val_adapter_fasta )
        ch_versions = ch_versions.mix( BLAST_MAKEBLASTDB.out.versions )

        //
        // ADAPTER SEARCH WITH BLASTN
        //
        // Convert reads to FASTA for BLASTN
        SEQKIT_FQ2FA ( ch_input_trim_branch.fastq )
        ch_versions = ch_versions.mix( SEQKIT_FQ2FA.out.versions )
        SAMTOOLS_FASTA ( ch_input_trim_branch.bam, false )

        fasta_for_blast = ch_input_trim_branch.fasta
            .mix( SEQKIT_FQ2FA.out.fasta )
            .mix( SAMTOOLS_FASTA.out.other )

        BLAST_BLASTN ( fasta_for_blast, BLAST_MAKEBLASTDB.out.db.collect(), [],[],[] )
        ch_versions = ch_versions.mix ( BLAST_BLASTN.out.versions )

        //
        // PROCESS BLAST OUTPUT WITH HIFITRIMMER PROCESSBLAST
        //
        // Prepare input for Hifitimmer processblast
        ch_input_processblast = BLAST_BLASTN.out.txt.combine( ch_adapter_yaml, by: 0 )
            .multiMap { meta, blastn, yaml ->
                blastn: [ meta, blastn ]
                yaml: [ meta, yaml ]
            }

        HIFITRIMMER_PROCESSBLAST ( ch_input_processblast.blastn, ch_input_processblast.yaml )

        hifitrimmer_summary = hifitrimmer_summary.mix ( HIFITRIMMER_PROCESSBLAST.out.summary )
        hifitrimmer_bed = hifitrimmer_bed.mix ( HIFITRIMMER_PROCESSBLAST.out.bed )

        //
        // FILTER READS WITH HIFITRIMMER FILTERBAM
        //
        // Convert FASTA and FASTQ to BAM for hifitrimmer filtering
        FQ2BAM ( ch_input_trim_branch.fasta.mix( ch_input_trim_branch.fastq ) )
        bam_for_hifitrimmer = FQ2BAM.out.bam.mix( ch_input_trim_branch.bam )

        ch_input_filterbam = bam_for_hifitrimmer.combine( HIFITRIMMER_PROCESSBLAST.out.bed, by: 0 )
        HIFITRIMMER_FILTERBAM ( ch_input_filterbam )

        trimmed_fastx =  trimmed_fastx.mix( HIFITRIMMER_FILTERBAM.out.filtered )

        if ( val_output_format == 'cram' ) {
            // convert INPUTs to CRAMs to export, need args `--output-format cram`
            FQ2CRAM_TRIM ( trimmed_fastx )
            trimmed_cram = trimmed_cram.mix( FQ2CRAM_TRIM.out.cram )
        }
    } else {
        ch_input_skip_trim = ch_input_to_trim
    }

    ch_input_skip_trim_branch = ch_input_skip_trim
        .branch { meta, reads -> def filename = reads.toString()
            bam: filename.endsWith('.bam')
                return [ meta, reads ]
            fastx: true
                return [ meta, reads ]
        }

    //
    // PROCESS UNTRIMMED READS FOR OUTPUT
    //
    if ( val_output_format == 'cram' ) {
        // convert INPUTs to CRAMs to export, need args `--output-format cram`
        FQ2CRAM_UNTRIM ( ch_input_skip_trim_branch.fastx )
        SAMTOOLS_VIEW ( ch_input_skip_trim_branch.bam.map{ meta, bam_file -> [ meta, bam_file, []]}, [[],[]], [], [] )
        untrimmed_cram = untrimmed_cram
            .mix( FQ2CRAM_UNTRIM.out.cram )
            .mix( SAMTOOLS_VIEW.out.cram )
    } else {
        // Convert reads to FASTQs to export
        // with --t arg for samtools fastq to copy RG, BC,  an QT tags to FASTQ header line
        SAMTOOLS_FASTQ ( ch_input_skip_trim_branch.bam, false )

        untrimmed_fastx = untrimmed_fastx
            .mix( SAMTOOLS_FASTQ.out.other )
            .mix( ch_input_skip_trim_branch.fastx )
    }


    emit:
    untrimmed_fastx     = untrimmed_fastx      // [meta, fastx] untrimmed reads in FASTA/FASTQ format, if trimming, only trimmed files
    untrimmed_cram      = untrimmed_cram       // [meta, cram] untrimmed reads in CRAM format, if trimming, only trimmed files
    trimmed_fastx       = trimmed_fastx        // [meta, fastx] preprocessed reads in FASTA/FASTQ format, if trimming, only trimmed files
    trimmed_cram        = trimmed_cram         // [meta, cram] preprocessed reads in CRAM format, if trimming, only trimmed files
    lima_report         = lima_reports         // [meta, report]
    lima_summary        = lima_summary         // [meta, summary]
    hifitrimmer_bed     = hifitrimmer_bed      // [meta, bed]
    hifitrimmer_summary = hifitrimmer_summary  // [meta, summary]
    pbmarkdup_stat      = pbmarkdup_stats      // [meta, log]
    versions            = ch_versions
}
