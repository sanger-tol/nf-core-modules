include { BLAST_BLASTN                         } from '../../../modules/nf-core/blast/blastn/main'
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
include { TABIX_BGZIP as BGZIP_BLASTN          } from '../../../modules/nf-core/tabix/bgzip/main'
include { UNTAR                                } from '../../../modules/nf-core/untar/main'

workflow PACBIO_PREPROCESS {

    take:
    ch_reads                // Channel [meta, input]: input reads in FASTA/FASTQ/BAM format, if trimming, only FASTQ/BAM
    ch_adapter_yaml         // Channel [meta, yaml]: yaml file for hifitrimmer adapter trimming
    adapter_db              // adapter database for blastn
    uli_primers             // Primer file for lima

    main:
    ch_versions = Channel.empty()

    lima_reports = Channel.empty()
    lima_summary = Channel.empty()
    pbmarkdup_stats = Channel.empty()
    if (uli_primers) {
        //
        // DEMULTIPLEX WITH LIMA
        //
        LIMA( ch_reads, uli_primers )
        ch_versions = ch_versions.mix( LIMA.out.versions )

        lima_reports = lima_reports.mix( LIMA.out.report )
        lima_summary = lima_summary.mix( LIMA.out.summary )

        //
        // MARK DUPLICATES WITH PBMARKDUP
        //
        ch_input_for_md = LIMA.out.bam
            .mix(LIMA.out.fastq)
            .mix(LIMA.out.fasta)
            .mix(LIMA.out.fastqgz)
            .mix( LIMA.out.fastagz )

        PBMARKDUP( ch_input_for_md )
        ch_versions = ch_versions.mix( PBMARKDUP.out.versions )
        pbmarkdup_stats = pbmarkdup_stats.mix( PBMARKDUP.out.log )

        ch_input_for_trimming = PBMARKDUP.out.markduped
    } else {
        ch_input_for_trimming = ch_reads
    }

    reads_to_filter = ch_input_for_trimming
        .branch { meta, reads ->
            def filename = reads.toString()
            fasta: filename =~ /\.(fasta|fa|fna)(\.gz)?$/
                return [ meta, reads ]
            fastq: filename =~ /\.(fastq|fq)(\.gz)?$/
                return [ meta, reads ]
            bam: filename.endsWith('.bam')
                return [ meta, reads ]
        }

    hifitrimmer_summary = Channel.empty()
    hifitrimmer_bed = Channel.empty()
    if (adapter_db && ch_adapter_yaml) {
        // UNTAR adapter database
        UNTAR( adapter_db )
        ch_versions = ch_versions.mix( UNTAR.out.versions )
        //
        // ADAPTER SEARCH WITH BLASTN
        //
        // Convert reads to FASTA for BLASTN
        SEQKIT_FQ2FA ( reads_to_filter.fastq )
        ch_versions = ch_versions.mix( SEQKIT_FQ2FA.out.versions )
        SAMTOOLS_FASTA ( reads_to_filter.bam, false )
        ch_versions = ch_versions.mix(SAMTOOLS_FASTA.out.versions)

        fasta_for_blast = reads_to_filter.fasta
            .mix( SEQKIT_FQ2FA.out.fasta )
            .mix( SAMTOOLS_FASTA.out.other )

        BLAST_BLASTN ( fasta_for_blast, UNTAR.out.untar.collect(), [],[],[] )
        BGZIP_BLASTN ( BLAST_BLASTN.out.txt )
        ch_versions = ch_versions.mix ( BLAST_BLASTN.out.versions )

        //
        // PROCESS BLAST OUTPUT WITH HIFITRIMMER PROCESSBLAST
        //
        // Prepare input for Hifitimmer processblast
        ch_input_processblast = BGZIP_BLASTN.out.output.combine( ch_adapter_yaml, by: 0 )
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
        FQ2BAM ( reads_to_filter.fasta.mix( reads_to_filter.fastq ) )
        ch_versions = ch_versions.mix ( FQ2BAM.out.versions )
        bam_for_hifitrimmer = FQ2BAM.out.bam.mix( reads_to_filter.bam )

        ch_input_filterbam = bam_for_hifitrimmer.combine( HIFITRIMMER_PROCESSBLAST.out.bed, by: 0 )
        HIFITRIMMER_FILTERBAM ( ch_input_filterbam )

        fastx =  HIFITRIMMER_FILTERBAM.out.filtered

        // convert INPUTs to CRAMs to export, need args `--output-format cram`
        FQ2CRAM_TRIM ( fastx )
        cram = FQ2CRAM_TRIM.out.cram
        ch_versions = ch_versions.mix ( FQ2CRAM_TRIM.out.versions )
    } else {
        // Convert reads to FASTQs to export
        // with --t arg for samtools fastq to copy RG, BC,  an QT tags to FASTQ header line
        SAMTOOLS_FASTQ ( reads_to_filter.bam, false )
        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

        fastx = SAMTOOLS_FASTQ.out.other
            .mix( reads_to_filter.fastq )
            .mix( reads_to_filter.fasta )

        // convert INPUTs to CRAMs to export, need args `--output-format cram`
        FQ2CRAM_UNTRIM ( reads_to_filter.fasta.mix( reads_to_filter.fastq ) )
        ch_versions = ch_versions.mix ( FQ2CRAM_UNTRIM.out.versions )

        SAMTOOLS_VIEW ( reads_to_filter.bam.map{ meta, bam_file -> [ meta, bam_file, []]}, [[],[]], [], [] )
        cram = FQ2CRAM_UNTRIM.out.cram.mix( SAMTOOLS_VIEW.out.cram )
    }

    emit:
    fastx               = fastx                // [meta, fastx] preprocessed reads in FASTA/FASTQ format, if trimming, only trimmed files
    cram                = cram                 // [meta, cram] preprocessed reads in CRAM format, if trimming, only trimmed files
    lima_report         = lima_reports         // [meta, report]
    lima_summary        = lima_summary         // [meta, summary]
    hifitrimmer_bed     = hifitrimmer_bed      // [meta, bed]
    hifitrimmer_summary = hifitrimmer_summary  // [meta, summary]
    pbmarkdup_stat      = pbmarkdup_stats      // [meta, log]
    versions            = ch_versions
}
