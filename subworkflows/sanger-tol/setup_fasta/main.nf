include { GUNZIP                        } from "../../../modules/nf-core/gunzip/main"
include { GAWK as GAWK_UPPER_SEQUENCE   } from '../../../modules/nf-core/gawk/main'
include { SAMTOOLS_FAIDX                } from "../../../modules/nf-core/samtools/faidx/main"

workflow SETUP_FASTA {
    take:
    ch_reference // channel.of( [meta], reference )
    val_sizes    // boolean

    main:

    //
    // LOGIC: SPLIT THE INPUT REFERENCES INTO THOSE THAT ARE ZIPPED OR UNZIPPED
    //        THIS ENSURES THAT ALL INPUT ARE UNZIPPED FOR DOWNSTREAM PROCESSING
    //
    ch_reference
            .branch { _meta, file ->
                zipped: file.name.endsWith('.gz')
                unzipped: !file.name.endsWith('.gz')
            }
            .set {ch_input}


    //
    // MODULE: UNZIP INPUTS IF NEEDED
    //
    GUNZIP (
        ch_input.zipped
    )


    //
    // LOGIC: MIX CHANELS WHICH MAY OR MAY NOT BE EMPTY INTO A SINGLE QUEUE CHANNEL
    //
    unzipped_input = channel.empty()

    unzipped_reference = unzipped_input
        .mix(ch_input.unzipped, GUNZIP.out.gunzip)


    //
    // MODULE: UPPERCASE THE REFERENCE SEQUENCE
    //         AWK = IF HEADER PASS ELSE CONVERT LINE TO UPPER CASE
    //
    ch_upper_sequence = channel.of('''\
        /^>/ {
            print; next
        } {
            print toupper(\$0)
        }'''.stripIndent())
        .collectFile(name: "uppercase_sequence.awk", cache: true)
        .collect()

    GAWK_UPPER_SEQUENCE(
        unzipped_reference,
        ch_upper_sequence,
        false,
    )


    //
    // MODULE: GENERATE INDEX OF REFERENCE FASTA
    //
    SAMTOOLS_FAIDX (
        GAWK_UPPER_SEQUENCE.out.output.map { meta, file -> [meta, file, []] },
        val_sizes
    )


    emit:
    reference = GAWK_UPPER_SEQUENCE.out.output
    fai       = SAMTOOLS_FAIDX.out.fai
    sizes     = SAMTOOLS_FAIDX.out.sizes
}
