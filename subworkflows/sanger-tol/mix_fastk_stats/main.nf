// Subworkflow imports
include { GENOME_STATISTICS } from '../../sanger-tol/genome_statistics/main.nf'

// Module imports
include { FASTK_FASTK       } from '../../../modules/nf-core/fastk/fastk/main'
include { FIND_CONCATENATE  } from '../../../modules/nf-core/find/concatenate/main'

workflow MIX_FASTK_STATS {
    take:
    ch_reference            // channel.of ( [meta], file )
    ch_haplotypes           // channel.of ( [meta], [haplo1, haplo2, haplo* ...])
    ch_longreads            // channel.of ( [meta], [files] )
    ch_kmer_dir             // channel.of ( [meta], file )
    ch_busco                // channel.of ( [meta], file )
    val_busco_lineage       // string
    val_existing_fastk      // boolean
    val_run_busco           // boolean



    main:
    //
    // LOGIC: BUILD PRIMARY VS HAPLOTYPE CROSS PER META
    //        IF HAP LIST IS > 1, MIX INTO PSEUDO_PRIMARY AND PSEUDO_HAPLOTYPES
    //        IF HAP LIST IS 1, USE PRIMARY AND HAPLOTYPE
    //        IF HAP LIST IS 0, USE ONLY PRIMARY + []
    //
    def hap_files = ch_haplotypes
        .map { _meta, file -> file }
        .collect()
        .map { files -> [files] ?: [] }

    def mixed_assemblies = ch_reference
        .combine(hap_files)
        .map { meta, primary, hap_list ->
            def hap_list_safe = hap_list ?: []
            def new_meta = meta + [real_primary: primary]
            if (hap_list_safe.size() > 1) {
                def all_files = [primary] + hap_list_safe
                all_files.collect { file ->
                    def others = all_files.findAll { all_files != file }
                    def meta_with_type = new_meta + [ori_id: meta.id, type: (file == primary ? "PRIMARY" : "HAP")]
                    [meta_with_type, file, others]
                }
            } else if (hap_list_safe.size() == 1) {
                def meta_with_type = new_meta + [type: "PRIMARY"]
                [[meta_with_type, primary, [hap_list[0]]]]
            } else {
                def meta_with_type = new_meta + [type: "PRIMARY"]
                [[meta_with_type, primary, []]]
            }
        }
        .flatMap { item -> item }
        .collect()
        .map { rows ->
            rows.collate(3)
                .findAll { all_rows -> all_rows.size() == 3 }
                .withIndex()
                .collect { row, idx ->
                    def (meta, file, collection) = row
                    [ meta + [hap_id: "${meta.id}_${meta.type}_${idx}"], file, collection ]
                }
        }
        .flatMap { row -> row }


    def concat_input = mixed_assemblies.map { meta, file, collection ->
        def new_meta = meta + [
            id: meta.hap_id + (collection.size() > 1 ? new File(file.toString()).name + "_pseudo_haplotype" : "_prim_hap"),
            pseudo_primary: file ]
        [ new_meta, collection ]
    }


    //
    // MODULE: MERGE MULTIPLE HAPLOTYPES INTO A PSEUDO "PRIMARY/HAPLOTYPE"
    //         TO GET PER HAP RESULTS IN GENOME_STATISTICS
    //
    FIND_CONCATENATE (
        concat_input
    )


    //
    // LOGIC: IF PRE-COMPUTED FASTK DATA EXISTS FINE FILES, ELSE RUN FASTK
    //
    if (val_existing_fastk) {
        //
        // LOGIC: FIND THE HIST AND KTAB FILES FOR FASTK DATABASE UPDATES
        //
        ch_fastk_data = ch_kmer_dir
            .map { meta, dir ->
                def files = new File(dir.toString()).listFiles() ?: []
                def ktab = files.findAll { all_files -> all_files.name ==~ /.*\.ktab(\..+)?$/ }.collect { ktabs -> file(ktabs) }
                def hist = files.findAll { all_files -> all_files.name ==~ /.*\.hist$/ }.collect { hists -> file(hists) }
                tuple(meta, hist, ktab)
            }

    } else {

        //
        // MODULE: RUN FASTK TO GENERATE KTAB AND HIST FILES
        //
        FASTK_FASTK (
            ch_longreads
        )

        ch_fastk_data = FASTK_FASTK.out.ktab
            .combine(FASTK_FASTK.out.hist, by: 0)
            .map { meta, ktab, hist ->
                tuple(meta, [hist], ktab)
            }
    }


    //
    // LOGIC: COMBINE CHANNELS SO THAT THEY ARE READY FOR GENOME_STATISTICS
    //
    gs_input = FIND_CONCATENATE.out.file_out
        .combine(ch_fastk_data)
        .multiMap { meta, haplotype, meta_fk, hist, ktabs ->
            fastk_channel: [
                meta, hist[0], ktabs, [], []
            ]
            refer_channel: [
                meta, meta.pseudo_primary, haplotype
            ]
            busco_channel: [
                meta, val_busco_lineage ?: "auto"
            ]
        }


    //
    // SUBWORKFLOW: RUN GENOME_STATISTICS FOR BUSCO, GFASTATS AND MERQURY STATS
    //              THIS WILL RUN THE PRIMARY or PSEUDOPRIMARY AS THE MAIN ASSEMBLY
    //
    GENOME_STATISTICS (
         gs_input.refer_channel,                                            // Needs to be in a [meta], asm1, asm2 - merqury can only take 2 assemblies so we trick it as above.
         gs_input.fastk_channel,                                            // format [ val(meta), fastk_hist_file, [fastk ktab files], [mat_fastk_ktabs ?: []], [pat_fastk_ktabs ?: []] ]
         gs_input.busco_channel,                                            // [ val(meta), string: busco_lineage ] - "lepidoptera_odb12" for specific lineage, "auto", "auto_euk", "auto_prok" // NO CSV OPTION, MEANS WE MIGHT NEED TO STOP THIS ONE AND RUN IT IN BTK ONLY?
         val_run_busco ? ch_busco.map{ _meta, file -> file } : []     // path to dir
    )


    // NOTE: REFORMAT THE CHANNEL INTO SOMETHING MORE SENSIBLE
    re_format_cat = FIND_CONCATENATE.out.file_out
        .map{ meta, haplotype ->
            def new_meta = meta.findAll { meta_array -> meta_array.key != 'pseudo_primary' }
            tuple( new_meta, meta.pseudo_primary, haplotype )
        }

    emit:
    gstats               = GENOME_STATISTICS.out.stats
    asmstats             = GENOME_STATISTICS.out.asmstats
    gfstats              = GENOME_STATISTICS.out.gfastats
    busco                = GENOME_STATISTICS.out.busco
    busco_batch_summary  = GENOME_STATISTICS.out.busco_batch_summary
    busco_summary_txt    = GENOME_STATISTICS.out.busco_summary_txt
    busco_summary_json   = GENOME_STATISTICS.out.busco_summary_json
    busco_log            = GENOME_STATISTICS.out.busco_log
    busco_directory      = GENOME_STATISTICS.out.busco_directory
    merqury              = GENOME_STATISTICS.out.merqury
    merqury_qv           = GENOME_STATISTICS.out.merqury_qv
    merqury_completeness = GENOME_STATISTICS.out.merqury_completeness
    merqury_phased_stats = GENOME_STATISTICS.out.merqury_phased_stats
    merqury_images       = GENOME_STATISTICS.out.merqury_images
    fastkdb              = ch_fastk_data
    cross_hap_collection = re_format_cat

}
