// SANGER_TOL SUBWORKFLOW IMPORTS
include { PAIRS_CREATE_CONTACT_MAPS } from "../../../subworkflows/sanger-tol/pairs_create_contact_maps/main.nf"

// SANGER_TOL MODULE IMPORTS
include { PRETEXTGRAPH              } from "../../../modules/sanger-tol/pretextgraph/main.nf"


workflow CONTACT_MAPS_INGEST_TRACKS {
    take:
    ch_pairs                    // [meta, pairs]
    ch_chrom_sizes              // [meta, sizes]
    ch_custom_order             // [meta, order]
    ch_gap_file                 // [meta, gap]
    ch_coverage                 // [meta, cov]
    ch_telo_file                // [meta, telo]
    ch_repeat_density           // [meta, repeat]
    val_build_pretext           // bool: build pretext map
    val_create_pretext_snapshot // bool: build snapshot
    val_build_cooler            // bool: build cooler
    val_build_juicer            // bool: build juicer
    val_cool_bin                // val: cooler cload parameter
    val_ingest_tracks           // bool: ingest accessory files into pretext

    main:

    //
    // SUBWORKFLOW: MAP THE PRETEXT FILE AND TAKE SNAPSHOT
    //
    PAIRS_CREATE_CONTACT_MAPS (
        ch_pairs,
        [[:],[]],
        ch_custom_order,
        true,
        true,
        false,
        false,
        []
    )


    //
    // MODULE: INGEST ACCESSORY FILES INTO PRETEXT BY DEFAULT
    //
    PRETEXT_INGEST (
        CREATE_MAPS.out.pretext.filter { !val_ingest_tracks },
        ch_gap_file,
        ch_coverage,
        ch_telo_file,
        ch_repeat_density
    )


    //
    // MODULE: CREATE SNAPSHOT PNG OF PRETEXTMAP
    //
    PRETEXT_SNAPSHOT (
        CREATE_MAPS.out.pretext.filter { !val_create_pretext_snapshot }
    )


    emit:
    pretext_map             = PAIRS_CREATE_CONTACT_MAPS.out.pretext
    pretext_map_with_tracks = PRETEXTGRAPH.out.pretext
    pretext_png             = PAIRS_CREATE_CONTACT_MAPS.out.pretext_png
    cooler                  = PAIRS_CREATE_CONTACT_MAPS.out.cool
    hic                     = PAIRS_CREATE_CONTACT_MAPS.out.hic
}
