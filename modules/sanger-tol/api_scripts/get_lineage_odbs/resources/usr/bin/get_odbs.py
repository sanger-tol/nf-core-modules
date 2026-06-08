#!/usr/bin/env python3

import argparse
import json
import os
import sys
from pathlib import Path

import requests
import urllib3

GOAT_API = "https://goat.genomehubs.org/api/v2/search?query=tax_lineage%28TAXID%29%20&result=taxon&includeEstimates=true&taxonomy=ncbi"
ENA_API = "https://www.ebi.ac.uk/ena/taxonomy/rest/v2/tax-id/TAXID?binomialOnly=false"


def parse_args(args=None):
    description = "Get ODB database value using NCBI API and BUSCO configuration file"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--taxid", help="TaxID for the species to retrieve ODBs for.")
    parser.add_argument(
        "--odb_version",
        help="Version of ODB to use",
        choices=["odb10", "odb12", "odb12.2", "all"],
        action="append",
    )
    parser.add_argument("--file_out", help="Output CSV file.")
    parser.add_argument("--odb_dir", help="Directory containing ODB files (Excluding the lineage subdirectory).")
    parser.add_argument(
        "--mapping_dir",
        help="Directory containing mapping files (BUSCO taxid to lineage mapping txt files).",
        type=str,
    )
    parser.add_argument(
        "--mode",
        help="Mode for ODB selection criteria",
        choices=["ancestral", "basal", "latest"],
        action="append",
    )
    parser.add_argument("--debug", help="Debug mode.", action="store_true", default=False)

    parser.add_argument(
        "--extra_lineages", help="Specified lineage ODBs to use. Minus the _odb{int}", nargs="+", default=None
    )
    parser.add_argument(
        "--basal_lineages",
        help="Basal lineages for BUSCO analysis",
        nargs="+",
        default=["bacteria", "archaea", "eukaryota"],
    )
    parser.add_argument("--version", action="version", version="%(prog)s 2.0")

    # Quick input validation
    output_args = parser.parse_args(args)

    if all(["ancestral", "latest"]) in output_args.mode:
        parser.error("Cannot use 'ancestral' and 'latest' at the same time.")

    if not output_args.mode and not output_args.specified_lineages:
        parser.error("Must have either modes or specified lineages.")

    return output_args


def make_dir(path):
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)


def get_http_request_json(url):
    """
    Robust wrapper for requests.get()
    """
    retry_strategy = urllib3.util.Retry(total=5, backoff_factor=0.1, status_forcelist=[429, 500, 502, 503, 504])
    adapter = requests.adapters.HTTPAdapter(max_retries=retry_strategy)
    session = requests.Session()
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    response = session.get(url)
    return response.json()


def goat_api_call(taxid: int, lineage_tax_ids_dict: dict[int, dict]):
    response = get_http_request_json(GOAT_API.replace("TAXID", str(taxid)))
    data = response["results"][0]["result"]["lineage"]
    query_tax_list = [int(taxid["taxon_id"]) for taxid in data]

    return {
        lineage_tax_ids_dict[taxon_id]["lineage"]: {
            "lineage": lineage_tax_ids_dict[taxon_id]["lineage"],
            "taxid": lineage_tax_ids_dict[taxon_id]["taxid"],
            "odb_string": lineage_tax_ids_dict[taxon_id]["odb_string"],
        }
        for taxon_id in query_tax_list
        if taxon_id in lineage_tax_ids_dict
    }


def ena_api_call(taxid: int, lineage_tax_ids_dict: dict[int, dict]):
    ena_response = get_http_request_json(ENA_API.replace("TAXID", str(taxid)))
    lineage_name = [lineage.lower() for lineage in ena_response["lineage"].split(";")]
    odb_lineage_names = [tax_lin["lineage"] for tax_lin in lineage_tax_ids_dict.values()]

    return {
        lineage_tax_ids_dict[int(y["taxid"])]["lineage"]: {
            "lineage": lineage_tax_ids_dict[int(y["taxid"])]["lineage"],
            "taxid": lineage_tax_ids_dict[int(y["taxid"])]["taxid"],
            "odb_string": lineage_tax_ids_dict[int(y["taxid"])]["odb_string"],
        }
        for i in set(lineage_name + odb_lineage_names)
        for _x, y in lineage_tax_ids_dict.items()
        if y["lineage"] == i
    }


def get_lineage_data(taxid: int, lineage_tax_ids_dict: dict[int, dict]) -> dict[int, dict[str, str]]:
    """
    Get lineage data for a given taxid
    """
    try:
        odb_arr: dict[int, dict] = goat_api_call(taxid, lineage_tax_ids_dict)
    except (requests.RequestException, KeyError, IndexError) as e:
        print(f"Warning: GOAT API call failed ({e}), falling back to ENA API", file=sys.stderr)
        try:
            odb_arr: dict[int, dict] = ena_api_call(taxid, lineage_tax_ids_dict)
        except (requests.RequestException, KeyError, IndexError) as ena_error:
            print(
                f"Error: Both GOAT and ENA servers could not be contacted. GOAT error: {e}. ENA error: {ena_error}",
                file=sys.stderr,
            )
            sys.exit(1)

    return odb_arr


def get_odb(mode, taxid, basal_lineages, extra_lineages, lineage_tax_ids_dict, odb_string, debug) -> dict:
    """
    Read the mapping between the BUSCO lineages and their taxon_id
    """
    odb_arr = get_lineage_data(taxid, lineage_tax_ids_dict)

    master_list = dict()

    if "ancestral" in mode:
        master_list.update(
            {
                f"{odb_arr[lineage]['taxid']}_{lineage}{odb_string}": {
                    "odb": lineage + odb_string,
                    "taxid": odb_arr[lineage]["taxid"],
                    "class": (
                        "latest"
                        if lineage == list(odb_arr)[0]
                        else "ancestral"
                        if "basal" in mode and lineage not in basal_lineages
                        else "basal"
                        if "basal" in mode and lineage in basal_lineages
                        else "ancestral"
                    ),
                }
                for lineage in odb_arr
            }
        )

    if "latest" in mode:
        first_key = list(odb_arr)[0]
        master_list[f"{odb_arr[first_key]['taxid']}_{odb_arr[first_key]['lineage']}{odb_string}"] = {
            "odb": odb_arr[first_key]["lineage"] + odb_string,
            "taxid": odb_arr[first_key]["taxid"],
            "class": "ancestral",
        }

    if "basal" in mode:
        for basal in basal_lineages:
            for x, y in lineage_tax_ids_dict.items():
                if basal == y["lineage"]:
                    if f"{y['taxid']}_{basal}{odb_string}" not in master_list.keys():
                        master_list[f"{y['taxid']}_{basal}{odb_string}"] = {
                            "odb": basal + odb_string,
                            "taxid": y["taxid"],
                            "class": "basal",
                        }

    if extra_lineages:
        for lineage in extra_lineages:
            for x, y in lineage_tax_ids_dict.items():
                if lineage == y["lineage"]:
                    if f"{y['taxid']}_{lineage}{odb_string}" not in master_list.keys():
                        master_list[f"{y['taxid']}_{lineage}{odb_string}"] = {
                            "odb": lineage + odb_string,
                            "taxid": y["taxid"],
                            "class": "extra",
                        }

    if debug:
        print(master_list)

    return master_list


def print_out(lineage_list, file_out, debug):
    """
    Print the lineage list to the output file.
    One line per lineage
    """
    with open(file_out, "w") as fout:
        for item_code, data in lineage_list.items():
            line = f"busco_lineage,{data['odb']},{data['class']}\n"
            if debug:
                print(line)

            fout.write(line)


def validate_lineage(lineage: dict, lineages_path: str):
    """
    Validate that the lineage exists in the lineage path.
    IF path is given, if not then we assume that the user want to run busco in ONLINE mode which means we can't validate local ODBs.
    """
    error_lineages = []
    for x, y in lineage.items():
        if lineages_path and not os.path.exists(lineages_path + "/lineages/" + y["odb"]):
            error_lineages.append(y["odb"])
        else:
            print(
                f"Skipping validation of {lineage}, odb_dir not provided (indicates your probably wanting busco to run in online mode)"
            )

    if len(error_lineages) > 0:
        raise FileNotFoundError(f"Lineages {error_lineages} not found in {lineages_path}")


def get_mapping_file(mapping_dir: str, odb_version: list, debug: bool):
    """
    Get the mapping file(s) for the ODB version.
    Returns a list of (mapping_file, odb_version_string) tuples.
    If odb_version is 'all', returns both odb10 and odb12.
    """
    mapping_files = []
    odb_version_list = ["odb10", "odb12", "odb12.2"] if "all" in odb_version else odb_version

    for file in Path(mapping_dir).glob("mapping_taxids-busco_dataset_name.eukaryota_*.txt"):
        for odb in odb_version_list:
            # There is no odb12.2 mapping file as of 5th June 2026
            odb_placeholder = "odb12" if odb == "odb12.2" else odb
            if f"_{odb_placeholder}" in str(file):
                (print("Found", file) if debug else None)
                mapping_files.append((str(file), f"_{odb}"))

    return mapping_files


def read_mapping_file(mapping_file: str, odb_string: str):
    """
    Read the mapping file and return a dictionary mapping lineage taxids
    to odb_dataset data.
    """
    lineage_tax_ids_dict = {}
    with open(mapping_file) as file_in:
        for line in file_in:
            arr = line.split()
            lineage_tax_ids_dict[int(arr[0])] = {
                "taxid": arr[0],
                "lineage": arr[1],
                "odb_string": odb_string,
            }
    return lineage_tax_ids_dict


def main(args=None):
    args = parse_args(args)

    mapping_files = get_mapping_file(args.mapping_dir, args.odb_version, args.debug)
    all_lineages = dict()

    for mapping_file, odb_version_string in mapping_files:
        mapping_data = read_mapping_file(mapping_file, odb_version_string)

        lineage_list = get_odb(
            args.mode,
            args.taxid,
            args.basal_lineages,
            args.extra_lineages,
            mapping_data,
            odb_version_string,
            args.debug,
        )

        all_lineages.update(sorted(lineage_list.items(), key=lambda x: x[1]["taxid"]))

    validate_lineage(all_lineages, args.odb_dir)
    make_dir(os.path.dirname(args.file_out))
    print_out(all_lineages, args.file_out, args.debug)


if __name__ == "__main__":
    main()
