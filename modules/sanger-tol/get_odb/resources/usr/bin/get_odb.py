#!/usr/bin/env python3

import argparse
import os
import json
import sys
import requests
import re


NCBI_TAXONOMY_API = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/%s"


def parse_args(args=None):
    Description = "Get ODB database value using NCBI API and BUSCO configuration file"

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("NCBI_SUMMARY_JSON", help="NCBI entry for this assembly for this assembly (in JSON).")
    parser.add_argument("LINEAGE_TAX_IDS", help="Mapping between BUSCO lineages and taxon IDs.")
    parser.add_argument("FILE_OUT", help="Output CSV file.")
    parser.add_argument("--all_ancestral_lineages", help="A boolean flag to include all ancestral lineages or just the best", action="store_true")
    parser.add_argument("--basal_lineages", help="Lineages considered ancestral for ODB selection.", nargs="+", default=["bacteria", "archaea", "eukaryota"])
    parser.add_argument("--version", action="version", version="%(prog)s 1.2")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)


def get_odb_version(file_name):
    match file_name:
        case _ if "odb10" in file_name:
            return "_odb10"
        case _ if "odb12" in file_name:
            return "_odb12"
        case _:
            sys.exit("Not a recognised ODB")


def get_odb(ncbi_summary, lineage_tax_ids, file_out, all_ancestral_lineages, basal_lineages):
    # Read the mapping between the BUSCO lineages and their taxon_id
    with open(lineage_tax_ids) as file_in:
        lineage_tax_ids_dict = {}
        for line in file_in:
            arr = line.split()
            lineage_tax_ids_dict[int(arr[0])] = arr[1]

    # Get the taxon_id of this species / assembly
    with open(ncbi_summary) as file_in:
        data = json.load(file_in)
    tax_id = data["reports"][0]["organism"]["tax_id"]

    # Using API, get the taxon_ids of all parents
    response = requests.get(NCBI_TAXONOMY_API % tax_id).json()
    ancestor_taxon_ids = response["taxonomy_nodes"][0]["taxonomy"]["lineage"]

    # Do the intersection to find the ancestors that have a BUSCO lineage
    odb_arr = [lineage_tax_ids_dict[taxon_id] for taxon_id in ancestor_taxon_ids if taxon_id in lineage_tax_ids_dict]

    # Get the ODB version from the file name
    odb_version = get_odb_version(lineage_tax_ids)

    # Include the basal lineages if needed, if all_ancestral_lineages is False, then we assume the user only wants lineage of best fit.
    odb_arr: list[str] = odb_arr + [ i for i in basal_lineages if not all_ancestral_lineages]

    # If all_ancestral_lineages is True, return all ancestral ODB lineages, otherwise return the closest ODB lineage
    if all_ancestral_lineages:
        odb_val =  [ i + odb_version for i in odb_arr]
    else:
        # The most recent [-1] OBD10/ODB12 lineage is selected
        odb_val = odb_arr[-1] + odb_version

    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)

    with open(file_out, "w") as fout:
        print("busco_lineage", odb_val, sep=",", file=fout)


def main(args=None):
    args = parse_args(args)
    get_odb(
        args.NCBI_SUMMARY_JSON,
        args.LINEAGE_TAX_IDS,
        args.FILE_OUT,
        args.all_ancestral_lineages,
        args.basal_lineages
    )


if __name__ == "__main__":
    sys.exit(main())
