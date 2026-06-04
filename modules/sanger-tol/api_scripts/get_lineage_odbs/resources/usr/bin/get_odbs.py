#!/usr/bin/env python3

import argparse
import json
import os
import sys

import requests

GOAT_API = "https://goat.genomehubs.org/api/v2/search?query=tax_lineage%28TAXID%29%20&result=taxon&includeEstimates=true&taxonomy=ncbi"
ENA_API = "https://www.ebi.ac.uk/ena/taxonomy/rest/v2/tax-id/TAXID?binomialOnly=false"


def parse_args(args=None):
    description = "Get ODB database value using NCBI API and BUSCO configuration file"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--taxid", help="TaxID for the species to retrieve ODBs for.")
    parser.add_argument("--odb_version", help="Version of ODB to use", choices=["all", "odb10", "odb12"])
    parser.add_argument("--file_out", help="Output CSV file.")
    parser.add_argument("--odb_dir", help="Directory containing ODB files (Excluding the lineage subdirectory).")
    parser.add_argument(
        "--mode",
        help="Mode for ODB selection criteria",
        choices=[
            "ancestral",
            "ancestral_and_basal",
            "best",
            "best_and_basal",
            "basal_only",
            "specified",
            "specified_and_basal",
        ],
        default="best",
    )
    parser.add_argument(
        "--specified_lineages", help="Specified lineage ODBs to use. Minus the _odb{int}", nargs="+", default=None
    )
    parser.add_argument(
        "--basal_lineages",
        help="Basal lineages for BUSCO analysis",
        nargs="+",
        default=["bacteria", "archaea", "eukaryota"],
    )
    parser.add_argument("--version", action="version", version="%(prog)s 1.3")

    # Quick input validation
    output_args = parser.parse_args(args)

    if output_args.mode in ["specified", "specified_and_basal"] and (output_args.specified_lineages is None):
        parser.error("--specified_lineages is required when mode is 'specified' or 'specified_and_basal'.")

    return output_args


def make_dir(path):
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)


def get_taxid(ncbi_summary: str) -> str:
    """
    Get taxon_id of the species from the NCBI summary JSON file.
    """
    with open(ncbi_summary) as file_in:
        data = json.load(file_in)
    return data["reports"][0]["organism"]["tax_id"]


def goat_api_call(taxid: int, lineage_tax_ids_dict: dict[int, str]):
    response = requests.get(GOAT_API.replace("TAXID", str(taxid))).json()
    data = response["results"][0]["result"]["lineage"]
    query_tax_list = [int(taxid["taxon_id"]) for taxid in data]
    odb_arr = [lineage_tax_ids_dict[taxon_id] for taxon_id in query_tax_list if taxon_id in lineage_tax_ids_dict]
    return odb_arr


def ena_api_call(taxid: int, lineage_tax_ids_dict: dict[int, str]):
    ena_reponse = requests.get(ENA_API.replace("TAXID", str(taxid))).json()
    lineage_name = [lineage.lower() for lineage in ena_reponse["lineage"].split(";")]
    odb_lineage_names = lineage_tax_ids_dict.values()
    odb_arr = [i for i in lineage_name if i in odb_lineage_names]
    return odb_arr


def get_lineage_data(taxid: int, lineage_tax_ids_dict: dict[int, str]) -> list[str]:
    """
    Get lineage data for a given taxid
    """
    try:
        odb_arr: list[str] = goat_api_call(taxid, lineage_tax_ids_dict)
    except (requests.RequestException, KeyError, IndexError) as e:
        print(f"Warning: GOAT API call failed ({e}), falling back to ENA API", file=sys.stderr)
        try:
            odb_arr: list[str] = ena_api_call(taxid, lineage_tax_ids_dict)
        except (requests.RequestException, KeyError, IndexError) as ena_error:
            print(
                f"Error: Both GOAT and ENA servers could not be contacted. GOAT error: {e}. ENA error: {ena_error}",
                file=sys.stderr,
            )
            sys.exit(1)

    return odb_arr


def get_odb(mode, taxid, odb_path, basal_lineages, mapping_file, odb_string):

    # Read the mapping between the BUSCO lineages and their taxon_id
    with open(mapping_file) as file_in:
        lineage_tax_ids_dict = {}
        for line in file_in:
            arr = line.split()
            lineage_tax_ids_dict[int(arr[0])] = arr[1]

    odb_arr = get_lineage_data(taxid, lineage_tax_ids_dict)

    if mode in ["ancestral", "ancestral_and_basal"]:
        # if user has requested ancestral lineages, return all of them up the tree
        odb_val: list[str] = [validate_lineage(lineage + odb_string, odb_path) for lineage in odb_arr]
    elif mode in ["best", "best_and_basal"]:
        # The most recent [-1] OBD10/ODB12 lineage is selected
        # In this case we only want the closest lineage, so we exclude the basal lineages
        # Force it into a list so that it isn't mangled by the interating print statement
        odb_val: list[str] = [validate_lineage(odb_arr[-1] + odb_string, odb_path)]

    if mode in ["ancestral_and_basal", "best_and_basal"]:
        # Add basals if the user wants them.
        odb_val: list[str] = odb_val + [lineage for lineage in basal_lineages if lineage not in odb_arr]

    return odb_val


def print_out(lineage_list, file_out):
    """
    Print the lineage list to the output file.
    One line per lineage
    """
    with open(file_out, "w") as fout:
        for lineage in lineage_list:
            print("busco_lineage", lineage, sep=",", file=fout)


def validate_lineage(lineage: str, lineages_path: str) -> str:
    """
    Validate that the lineage exists in the lineage path.
    IF path is given, if not then we assume that the user want to run busco in ONLINE mode which means we can't validate local ODBs.
    """
    if lineages_path:
        if not os.path.exists(lineages_path + "/lineages/" + lineage):
            raise FileNotFoundError(f"Lineage {lineage} not found in {lineages_path}")
    else:
        print(
            "Skipping validation of ODB files, as odb_dir not provided (indicates your probably wanting busco to run in online mode"
        )
    return lineage


def get_specific_odbs(
    lineage_path,
    specified_lineages,
    odb_version,
    basal_lineages=[],
):
    """
    User provided lineages will have the ODB
    """
    specified_odbs = [validate_lineage(i + odb_version, lineage_path) for i in specified_lineages]

    return specified_odbs + basal_lineages


def get_mapping_file(odb_version: str):
    """
    Get the mapping file(s) for the ODB version.
    Returns a list of (mapping_file, odb_version_string) tuples.
    If odb_version is 'all', returns both odb10 and odb12.
    """
    current_working_dir = os.path.dirname(os.path.realpath(__file__))
    cwd = current_working_dir + "/../assets/"

    mapping: dict[str, dict[str, str]] = {
        "odb10": {"file": "mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt", "version": "_odb10"},
        "odb12": {"file": "mapping_taxids-busco_dataset_name.eukaryota_odb12.2025-01-15.txt", "version": "_odb12"},
    }

    if odb_version == "all":
        return [(f"{cwd}/{mapping[v]['file']}", mapping[v]["version"]) for v in ["odb10", "odb12"]]
    else:
        odb_version_string = mapping[odb_version]["version"]
        mapping_file = f"{cwd}/{mapping[odb_version]['file']}"
        return [(mapping_file, odb_version_string)]


def main(args=None):
    args = parse_args(args)

    mapping_files = get_mapping_file(args.odb_version)
    all_lineages = []

    for mapping_file, odb_version_string in mapping_files:
        if "basal" in args.mode:
            new_basal_lineages: list[str] = [
                validate_lineage(i + odb_version_string, args.odb_dir) for i in args.basal_lineages
            ]
        else:
            new_basal_lineages = []

        if args.specified_lineages is not None and args.mode in ["specified", "specified_and_basal"]:
            lineage_list = get_specific_odbs(
                args.odb_dir,
                args.specified_lineages,
                odb_version_string,
                new_basal_lineages if args.mode == "specified_and_basal" else [],
            )

        elif args.mode in ["ancestral", "ancestral_and_basal", "best", "best_and_basal"]:
            lineage_list = get_odb(
                args.mode,
                args.taxid,
                args.odb_dir,
                new_basal_lineages if "basal" in args.mode else [],
                mapping_file,
                odb_version_string,
            )

        elif args.mode == "basal_only":
            lineage_list = get_specific_odbs(
                args.odb_dir,
                args.basal_lineages,
                odb_version_string,
                None,
            )

        all_lineages.extend(lineage_list)

    # Sort the list by the odb version so it looks nicer
    all_lineages = sorted(list(set(all_lineages)))

    make_dir(os.path.dirname(args.file_out))
    print_out(all_lineages, args.file_out)


if __name__ == "__main__":
    main()
