#!/usr/bin/env python3

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from string import Template

import requests
import urllib3

GOAT_API = Template(
    "https://goat.genomehubs.org/api/v2/search?query=tax_lineage%28${taxid}%29&result=taxon&includeEstimates=true&taxonomy=ncbi"
)
ENA_API = Template("https://www.ebi.ac.uk/ena/taxonomy/rest/v2/tax-id/${taxid}?binomialOnly=false")


@dataclass
class BuscoLineage:
    taxid: int
    lineage: str
    odb_string: str


@dataclass
class BuscoSelection:
    odb_string: str
    taxid: int
    classification: str


def parse_args(args=None):
    description = "Get ODB database value using NCBI API and BUSCO configuration file"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--taxid", help="TaxID for the species to retrieve ODBs for.", type=int, required=True)
    parser.add_argument(
        "--odb_version",
        help="Version of ODB to use",
        choices=["odb10", "odb12", "odb12.2", "all"],
        action="append",
    )
    parser.add_argument("--file_out", help="Output CSV file.", type=str, required=True)
    parser.add_argument("--odb_dir", help="Directory containing ODB files (Excluding the lineage subdirectory).")
    parser.add_argument(
        "--mapping_dir",
        help="Directory containing mapping files (BUSCO taxid to lineage mapping txt files).",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--mode",
        help="Mode for ODB selection criteria",
        choices=["ancestral", "basal", "latest", "none"],
        action="append",
        default=["none"],
    )
    parser.add_argument("--debug", help="Debug mode.", action="store_true", default=False)
    parser.add_argument(
        "--print_output",
        help="Print output to stdout. Doesn't print the full debug log.",
        action="store_true",
        default=False,
    )
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

    if not output_args.odb_version:
        parser.error("At least one --odb_version must be specified.")

    if all(["ancestral", "latest"]) in output_args.mode:
        parser.error("Cannot use 'ancestral' and 'latest' at the same time.")

    if not output_args.mode and not output_args.specified_lineages:
        parser.error("Must have either modes or specified lineages.")

    return output_args


def make_dir(path):
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)


def get_http_request_json(url: str):
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


def goat_api_call(taxid: int, lineage_tax_ids_dict: dict[int, BuscoLineage]) -> dict[str, BuscoLineage]:
    response = get_http_request_json(GOAT_API.substitute(taxid=taxid))
    data = response["results"][0]["result"]["lineage"]
    query_tax_list = [int(taxid["taxon_id"]) for taxid in data]

    return {
        lineage_tax_ids_dict[taxon_id].lineage: lineage_tax_ids_dict[taxon_id]
        for taxon_id in query_tax_list
        if taxon_id in lineage_tax_ids_dict
    }


def ena_api_call(taxid: int, lineage_tax_ids_dict: dict[int, BuscoLineage]) -> dict[str, BuscoLineage]:
    ena_response = get_http_request_json(ENA_API.substitute(taxid=taxid))
    lineage_name = [lineage.lower() for lineage in ena_response["lineage"].split(";")]
    odb_lineage_names = [tax_lin.lineage for tax_lin in lineage_tax_ids_dict.values()]

    return {
        lineage_tax_ids_dict[y.taxid].lineage: lineage_tax_ids_dict[y.taxid]
        for i in set(lineage_name + odb_lineage_names)
        for _x, y in lineage_tax_ids_dict.items()
        if y.lineage == i
    }


def get_lineage_data(taxid: int, lineage_tax_ids_dict: dict[int, BuscoLineage]) -> dict[str, BuscoLineage]:
    """
    Get lineage data for a given taxid
    """
    try:
        return goat_api_call(taxid, lineage_tax_ids_dict)
    except (requests.RequestException, KeyError, IndexError) as e:
        print(f"Warning: GOAT API call failed ({e}), falling back to ENA API", file=sys.stderr)
        try:
            return ena_api_call(taxid, lineage_tax_ids_dict)
        except (requests.RequestException, KeyError, IndexError) as ena_error:
            print(
                f"Error: Both GOAT and ENA servers could not be contacted. GOAT error: {e}. ENA error: {ena_error}",
                file=sys.stderr,
            )
            sys.exit(1)


def get_odb(
    mode: list[str],
    taxid: int,
    basal_lineages: list[str],
    extra_lineages: list[str],
    lineage_tax_ids_dict: dict[int, BuscoLineage],
    odb_string: str,
    debug: bool,
) -> dict[str, BuscoSelection]:
    """
    Read the mapping between the BUSCO lineages and their taxon_id
    """
    odb_dict: dict[str, BuscoLineage] = get_lineage_data(taxid, lineage_tax_ids_dict)

    master_list: dict[str, BuscoSelection] = dict()

    if "ancestral" in mode:
        master_list.update(
            {
                f"{odb_dict[lineage].taxid}_{lineage}{odb_string}": BuscoSelection(
                    odb_string=lineage + odb_string,
                    taxid=odb_dict[lineage].taxid,
                    classification=(
                        "latest"
                        if lineage == list(odb_dict)[0]
                        else "ancestral"
                        if "basal" in mode and lineage not in basal_lineages
                        else "basal"
                        if "basal" in mode and lineage in basal_lineages
                        else "ancestral"
                    ),
                )
                for lineage in odb_dict
            }
        )

    if "latest" in mode:
        first_key = list(odb_dict)[0]
        master_list[f"{odb_dict[first_key].taxid}_{odb_dict[first_key].lineage}{odb_string}"] = BuscoSelection(
            odb_string=odb_dict[first_key].lineage + odb_string,
            taxid=odb_dict[first_key].taxid,
            classification="ancestral",
        )

    if "basal" in mode:
        for basal in basal_lineages:
            for x, y in lineage_tax_ids_dict.items():
                if basal == y.lineage:
                    if f"{y.taxid}_{basal}{odb_string}" not in master_list.keys():
                        master_list[f"{y.taxid}_{basal}{odb_string}"] = BuscoSelection(
                            odb_string=basal + odb_string,
                            taxid=y.taxid,
                            classification="basal",
                        )

    if extra_lineages:
        for lineage in extra_lineages:
            for x, y in lineage_tax_ids_dict.items():
                if lineage == y.lineage:
                    if f"{y.taxid}_{lineage}{odb_string}" not in master_list.keys():
                        master_list[f"{y.taxid}_{lineage}{odb_string}"] = BuscoSelection(
                            odb_string=lineage + odb_string,
                            taxid=y.taxid,
                            classification="extra",
                        )

    if debug:
        print(master_list)

    return master_list


def print_out(lineage_list: dict[str, BuscoSelection], file_out: str, debug: bool):
    """
    Print the lineage list to the output file.
    One line per lineage
    """
    with open(file_out, "w") as fout:
        for item_code, data in lineage_list.items():
            line = f"{data.odb_string},{data.classification}"
            if debug:
                print(line)

            fout.write(f"{line}\n")


def validate_lineage(lineage: dict[str, BuscoSelection], lineages_path: str):
    """
    Validate that the lineage exists in the lineage path.
    IF path is given, if not then we assume that the user want to run busco in ONLINE mode which means we can't validate local ODBs.
    """
    error_lineages = []
    for data in lineage.values():
        if not os.path.exists(os.path.join(lineages_path, "lineages", data.odb_string)):
            error_lineages.append(data.odb_string)

    if len(error_lineages) > 0:
        raise FileNotFoundError(f"Lineages {error_lineages} not found in {lineages_path}")


def get_mapping_file(mapping_dir: str, odb_version: list, debug: bool) -> list[tuple[str, str]]:
    """
    Get the mapping file(s) for the ODB version.
    Returns a list of (mapping_file, odb_version_string) tuples.
    If odb_version is 'all', returns both odb10 and odb12.
    """
    mapping_files = []
    odb_version_list = ["odb10", "odb12", "odb12.2"] if "all" in odb_version else odb_version

    files = [str(file) for file in Path(mapping_dir).glob("*.txt")]
    for odb in odb_version_list:
        # There is no odb12.2 mapping file as of 5th June 2026
        # If a file is introduced we may have to change this logic to avoid matching
        # "mapping_taxids-busco_dataset_name.eukaryota_odb12.2025-01-15.txt" by accident
        # when looking for odb12.2 files
        odb_placeholder = "odb12" if odb == "odb12.2" else odb
        for file in files:
            if odb_placeholder in file:
                if debug:
                    print("Found", file)
                mapping_files.append((file, f"_{odb}"))

    if len(mapping_files) == 0:
        raise FileNotFoundError(f"No mapping files found in {mapping_dir} for odb version(s) {odb_version}")

    return mapping_files


def read_mapping_file(mapping_file: str, odb_string: str) -> dict[int, BuscoLineage]:
    """
    Read the mapping file and return a dictionary mapping lineage taxids
    to odb_dataset data.
    """
    lineage_tax_ids_dict = {}
    with open(mapping_file) as file_in:
        for line in file_in:
            arr = line.split()
            lineage_tax_ids_dict[int(arr[0])] = BuscoLineage(
                taxid=int(arr[0]),
                lineage=arr[1],
                odb_string=odb_string,
            )
    return lineage_tax_ids_dict


def main(args=None):
    args = parse_args(args)

    mapping_files = get_mapping_file(args.mapping_dir, args.odb_version, args.debug)
    all_lineages: dict[str, BuscoSelection] = dict()

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

        all_lineages.update(lineage_list.items())

    if args.odb_dir:
        validate_lineage(all_lineages, args.odb_dir)
    else:
        print("Skipping validation, odb_dir not provided (indicates your probably wanting busco to run in online mode)")

    make_dir(os.path.dirname(args.file_out))
    print_out(all_lineages, args.file_out, args.debug or args.print_output)


if __name__ == "__main__":
    main()
