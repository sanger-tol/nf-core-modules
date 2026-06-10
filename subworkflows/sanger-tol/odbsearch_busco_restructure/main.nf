include { API_SCRIPTS_GET_LINEAGE_ODBS  } from "../../../modules/sanger-tol/api_scripts/get_lineage_odbs"
include { BUSCO_BUSCO                   } from '../../../modules/sanger-tol/busco/busco/main'
include { RESTRUCTUREBUSCODIR           } from '../../../modules/sanger-tol/restructurebuscodir/main'

workflow ODBSEARCH_BUSCO_RESTRUCTURE {
    take:
    ch_reference                // tuple([meta], reference)
    val_odb_directory           // val(path to lineages folder)
    val_mapping_directory       // val(path to mapping folder)
    val_taxid                   // val(9606)
    val_specified_lineages      // val("mammalia")
    val_output_dir              // val(output directory)
    val_restructure_busco_dir   // val(boolean)

    main:

    //
    // MODULE: GET LIKELY ODB CANDIDATES FROM GOAT/ENA USING TAXID
    //
    API_SCRIPTS_GET_LINEAGE_ODBS (
        ch_reference,
        val_odb_directory,
        val_mapping_directory,
        val_taxid,
        val_specified_lineages
    )


    //
    // LOGIC: TAKE CSV AND RETURN FLAT LIST OF ODBS
    //
    ch_busco_input = API_SCRIPTS_GET_LINEAGE_ODBS.out.csv
        .map { meta, csv ->
            def odbs = csv
                .splitCsv( header: false )
                .collect { row -> [ row[1], row[2] ] }
            [ meta, odbs ]
        }
        .transpose() // Convert to [meta, odbs] pairs from [meta, [odb_list]]
        .map { meta, odb -> [ meta.id, meta, odb ] }
        .combine(
            ch_reference.map { meta, ref -> [ meta.id, meta, ref ] }, // Normalise the fasta so we can combine easier
            by: 0
        )
        .map { id, meta, odb, ref_meta, ref -> [ ref_meta, odb, ref ] }
        .combine( val_taxid )
        .combine( val_output_dir )
        .map { meta, odb, ref, tax_id, outdir_location ->
            def new_meta = meta + [ lineage: odb[0], lineage_rating: odb[1], taxid: tax_id, outdir: outdir_location, genome_size: ref.size() ]
            [ new_meta, ref ]
        }


    //
    // MODULE: Run BUSCO search
    //
    BUSCO_BUSCO(
        ch_busco_input,
        'genome',
        ch_busco_input.map { meta, _fasta -> meta.lineage },
        val_odb_directory,
        [],
        []
    )


    //
    // MODULE: Tidy up the BUSCO output directories before publication
    //
    RESTRUCTUREBUSCODIR(
        BUSCO_BUSCO.out.restructure_busco_output
    )


    emit:
    odb_csv             = API_SCRIPTS_GET_LINEAGE_ODBS.out.csv
    busco_full_table    = BUSCO_BUSCO.out.full_table
    busco_output        = BUSCO_BUSCO.out.busco_dir
    restructured_output = RESTRUCTUREBUSCODIR.out.clean_busco_dir
}
