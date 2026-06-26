include { APISCRIPTS_GETLINEAGEODBS     } from '../../../modules/sanger-tol/apiscripts/getlineageodbs'
include { BUSCO_BUSCO                   } from '../../../modules/sanger-tol/busco/busco/main'
include { RESTRUCTUREBUSCODIR           } from '../../../modules/sanger-tol/restructurebuscodir/main'

workflow ODBSEARCH_BUSCO_RESTRUCTURE {
    take:
    ch_reference                // tuple([meta], reference)
    val_odb_directory           // val(path to directory containing `lineages` folder)
    val_mapping_directory       // val(path to busco_odb_mapping folder, shipped with APISCRIPTS_GETLINEAGEODBS)
    ch_taxid                    // tuple([meta], val(9606))
    ch_specified_lineages       // tuple([meta], val("mammalia"))
    val_restructure_busco_dir   // val(boolean)

    main:

    //
    // MODULE: GET LIKELY ODB CANDIDATES FROM GOAT/ENA USING TAXID
    //
    APISCRIPTS_GETLINEAGEODBS (
        ch_reference,
        val_odb_directory,
        val_mapping_directory,
        ch_taxid.map{ meta, taxid -> taxid },
        ch_specified_lineages.map{ meta, taxid -> taxid }
    )


    //
    // LOGIC: TAKE CSV AND RETURN FLAT LIST OF ODBS
    //
    ch_busco_input = APISCRIPTS_GETLINEAGEODBS.out.csv
        .map { meta, csv ->
            def odbs = csv
                .splitCsv( header: false )
                .collect { row -> [ row[0], row[1] ] }
            [ meta, odbs ]
        }
        .transpose() // Convert to [meta, odbs] pairs from [meta, [odb_list]]
        .map { meta, odb -> [ meta.id, meta, odb ] }
        .combine(
            ch_reference.map { meta, ref -> [ meta.id, meta, ref ] }, // Normalise the fasta so we can combine easier
            by: 0,
        )
        .map { id, meta, odb, ref_meta, ref -> [ ref_meta, odb, ref ] }
        .unique { meta, odb, ref ->
            [ meta, odb ]
        } // Make unique by meta.id and odb[0] to avoid duplicate entries caused by multiple entried in the input samplesheet
        .multiMap { meta, odb, ref ->
            def new_meta = meta + [ lineage: odb[0], lineage_rating: odb[1] ]
            reference: [ new_meta, ref ]
            busco_mode: 'genome'
            lineage: new_meta.lineage
        }


    //
    // MODULE: Run BUSCO search
    //
    BUSCO_BUSCO(
        ch_busco_input.reference,
        ch_busco_input.busco_mode,
        ch_busco_input.lineage,
        val_odb_directory,
        [],
        []
    )


    //
    // MODULE: Tidy up the BUSCO output directories before publication
    //
    RESTRUCTUREBUSCODIR(
        BUSCO_BUSCO.out.restructure_busco_output.filter { val_restructure_busco_dir }
    )


    emit:
    odb_csv             = APISCRIPTS_GETLINEAGEODBS.out.csv
    busco_full_table    = BUSCO_BUSCO.out.full_table
    busco_output        = BUSCO_BUSCO.out.busco_dir
    restructured_output = RESTRUCTUREBUSCODIR.out.clean_busco_dir
}
