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


# A single entry from a mapping file
@dataclass(frozen=True)
class BuscoLineage:
    taxid: int
    lineage: str


# All entries from a mapping file, with helper functions to index them
@dataclass
class BuscoDatabase:
    lineages: list[BuscoLineage]

    def by_taxid(self) -> dict[int, BuscoLineage]:
        if not hasattr(self, "_lineage_tax_ids_dict"):
            self._lineage_tax_ids_dict = {lineage.taxid: lineage for lineage in self.lineages}
        return self._lineage_tax_ids_dict

    def by_lineage(self) -> dict[str, BuscoLineage]:
        if not hasattr(self, "_lineage_names_dict"):
            self._lineage_names_dict = {lineage.lineage: lineage for lineage in self.lineages}
        return self._lineage_names_dict

    def validate_and_get(self, name: str) -> BuscoLineage:
        try:
            return self.by_lineage()[name]
        except KeyError:
            raise ValueError(f"Lineage {name} not found in mapping file.")


# A class to hold the selected ODB lineages and their classifications (ancestral, basal, latest, extra)
class BuscoSelection(dict[str, str]):
    # Convenient method to add a lineage to the selection
    def add_lineage(self, lineage: BuscoLineage, classification: str, odb_string: str):
        self.setdefault(f"{lineage.lineage}{odb_string}", classification)


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
        choices=["ancestral", "basal", "latest"],
        action="append",
        default=[],
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

    if "ancestral" in output_args.mode and "latest" in output_args.mode:
        parser.error("Cannot use 'ancestral' and 'latest' at the same time.")

    if not output_args.mode and not output_args.extra_lineages:
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


def goat_api_call(taxid: int, lineage_tax_ids_dict: dict[int, BuscoLineage]) -> list[BuscoLineage]:
    response = get_http_request_json(GOAT_API.substitute(taxid=taxid))
    data = response["results"][0]["result"]["lineage"]
    query_tax_list = [int(taxid["taxon_id"]) for taxid in data]

    return [lineage_tax_ids_dict[taxon_id] for taxon_id in query_tax_list if taxon_id in lineage_tax_ids_dict]


def ena_api_call(taxid: int, lineage_names_dict: dict[str, BuscoLineage]) -> list[BuscoLineage]:
    ena_response = get_http_request_json(ENA_API.substitute(taxid=taxid))
    lineage_names = [lineage.strip().lower() for lineage in ena_response["lineage"].split(";")]

    return [lineage_names_dict[name] for name in reversed(lineage_names) if name in lineage_names_dict]


def get_lineage_data(taxid: int, busco_db: BuscoDatabase) -> list[BuscoLineage]:
    """
    Get lineage data for a given taxid
    """
    try:
        return goat_api_call(taxid, busco_db.by_taxid())
    except (requests.RequestException, KeyError, IndexError) as e:
        print(f"Warning: GOAT API call failed ({e}), falling back to ENA API", file=sys.stderr)
        try:
            return ena_api_call(taxid, busco_db.by_lineage())
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
    busco_db: BuscoDatabase,
    odb_string: str,
) -> BuscoSelection:
    """
    Read the mapping between the BUSCO lineages and their taxon_id
    """
    if "ancestral" in mode or "latest" in mode:
        ancestral_lineages = get_lineage_data(taxid, busco_db)
    else:
        ancestral_lineages = []

    master_list = BuscoSelection()

    if "ancestral" in mode:
        for i, lineage in enumerate(ancestral_lineages):
            if i == 0:
                classification = "latest"
            elif "basal" in mode and lineage.lineage in basal_lineages:
                classification = "basal"
            else:
                classification = "ancestral"
            master_list.add_lineage(lineage, classification, odb_string)

    if "latest" in mode:
        master_list.add_lineage(ancestral_lineages[0], "latest", odb_string)

    if "basal" in mode:
        for basal in basal_lineages:
            lin = busco_db.validate_and_get(basal)
            master_list.add_lineage(lin, "basal", odb_string)

    if extra_lineages:
        for lineage in extra_lineages:
            lin = busco_db.validate_and_get(lineage)
            master_list.add_lineage(lin, "extra", odb_string)

    return master_list


def print_out(lineage_list: BuscoSelection, file_out: str, debug: bool):
    """
    Print the lineage list to the output file.
    One line per lineage
    """
    with open(file_out, "w") as fout:
        for odb_string, classification in lineage_list.items():
            line = f"{odb_string},{classification}"
            if debug:
                print(line)

            fout.write(f"{line}\n")


def check_offline_availability(selected_buscos: BuscoSelection, lineages_path: str):
    """
    Validate that the lineage exists in the lineage path.
    IF path is given, if not then we assume that the user want to run busco in ONLINE mode which means we can't validate local ODBs.
    """
    error_lineages = []
    for odb_string in selected_buscos:
        if not os.path.exists(os.path.join(lineages_path, "lineages", odb_string)):
            error_lineages.append(odb_string)

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


def read_mapping_file(mapping_file: str) -> BuscoDatabase:
    """
    Read the mapping file and return the database object containing
    the list of taxids and lineage names.
    """
    lineages = []
    with open(mapping_file) as file_in:
        for line in file_in:
            arr = line.split()
            lineages.append(
                BuscoLineage(
                    taxid=int(arr[0]),
                    lineage=arr[1],
                )
            )
    return BuscoDatabase(lineages=lineages)


def main(args=None):
    args = parse_args(args)

    mapping_files = get_mapping_file(args.mapping_dir, args.odb_version, args.debug)
    all_lineages = BuscoSelection()

    for mapping_file, odb_version_string in mapping_files:
        mapping_data = read_mapping_file(mapping_file)

        lineage_list = get_odb(
            args.mode,
            args.taxid,
            args.basal_lineages,
            args.extra_lineages,
            mapping_data,
            odb_version_string,
        )

        if args.debug:
            print(odb_version_string, lineage_list)

        all_lineages.update(lineage_list)

    if args.odb_dir:
        check_offline_availability(all_lineages, args.odb_dir)
    else:
        print("Skipping validation, odb_dir not provided (indicates your probably wanting busco to run in online mode)")

    make_dir(os.path.dirname(args.file_out))
    print_out(all_lineages, args.file_out, args.debug or args.print_output)


if __name__ == "__main__":
    main()
