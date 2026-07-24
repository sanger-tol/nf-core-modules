#!/usr/bin/env python3
"""
FASTA sanitisation script
---------------------------

Processes FAI file to:
1. Read FAI and create dict with sequence name, length, status and new_name
2. Mark sequences as TRIM or KEEP based on sequence length threshold
3. Sanitise headers and add new names to dict
4. Output mapping file (old_name -> new_name) for KEEP sequences
5. Output a samtools header subsetting file.
5. Output full stats as JSON

Downstream is:
    AWK to replace old headers with new ones
    SAMTOOLS faidx to create a subset and indexed fasta/fai
    AWK to correct the sequence if needed

Written by DLBPointon based on:
    - filter_fasta_by_length
    - sanitise_fasta_header
"""

import argparse
import json
import re
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sanitise FASTA headers and filter by length using FAI file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("fai_file", type=str, help="Path to FAI file")
    parser.add_argument(
        "--max_sequence_length",
        type=int,
        default=1000000,
        help="A sequence cut off value, sequences LONGER than this are marked TRIM, others KEEP",
    )
    parser.add_argument(
        "--max_header_length",
        type=int,
        default=50,
        help="Maximum length of FASTA header",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="filter_fasta",
        help="Output prefix for naming files",
    )
    parser.add_argument("-v", "--version", action="version", version="1.0.0")
    return parser.parse_args()


def read_fai(fai_path: str) -> dict[str, dict[str, int | str]]:
    """
    Read FAI file and return dict of data
    """
    fai_dict = {}

    if not Path(fai_path).exists():
        raise FileNotFoundError(f"FAI file not found: {fai_path}")

    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split()
            if not (len(parts) == 2 or len(parts) == 5):
                sys.stderr.write(f"Are you sure this is a valid FAI/Sizes file? Offender: {fai_path}\n")
                sys.exit(127)

            name = parts[0]
            length = int(parts[1])

            if name in fai_dict:
                sys.stderr.write(f"Warning: Duplicate name {name} found, skipping\n")
                sys.exit(127)
            else:
                fai_dict[name] = {
                    "length": length,
                    "status": "N",
                    "new_name": "N",
                }

    return fai_dict


def assign_status(fai_dict: dict[str, dict[str, int | str]], threshold: float) -> dict[str, dict[str, int | str]]:
    """
    Assign TRIM or KEEP status based on length threshold
    """
    for name, data in fai_dict.items():
        if int(data["length"]) > int(threshold):
            data["status"] = "TRIM"
        else:
            data["status"] = "KEEP"

    return fai_dict


def sanitise_headers(fai_dict: dict, max_header_length: int) -> dict[str, dict[str, int | str]]:
    """
    Sanitise headers and add new_name to dict
    """
    counter = 0
    for name, data in fai_dict.items():
        counter += 1
        # Replace problematic characters with underscores
        sanitised = re.sub(r"[,;\s|:]", "_", name)

        # Truncate to max_header_length - 5 and append a counter
        # to avoid exceeding max length while keeping room for the counter
        # some headers are very long and identicle until past the max length
        if len(sanitised) > max_header_length:
            sanitised = sanitised[: max_header_length - 5] + f"_{counter:04d}"

        data["new_name"] = sanitised

    return fai_dict


def write_mapping_file(fai_dict: dict, output_path: str) -> None:
    """
    Write mapping file (old_name\\tnew_name) for KEEP sequences only
    """
    with open(output_path, "w") as f:
        for name, data in fai_dict.items():
            if data["status"] == "KEEP":
                print(f"{name}\t{data['new_name']}", file=f)


def write_new_scaffolds(fai_dict: dict, output_path: str) -> None:
    """
    Write new scaffold names to TXT file for subsetting a fasta by
    SAMTOOLS faidx -r
    """
    with open(output_path, "w") as f:
        for name, data in fai_dict.items():
            if data["status"] == "KEEP":
                print(data["new_name"], file=f)


def write_json_stats(fai_dict: dict, max_sequence_length: int, max_header_length: int, output_path: str) -> None:
    """
    Write full stats JSON file
    """
    stats = {
        "total_sequences": len(fai_dict),
        "keep_count": sum(1 for d in fai_dict.values() if d["status"] == "KEEP"),
        "trim_count": sum(1 for d in fai_dict.values() if d["status"] == "TRIM"),
        "max_sequence_length": max_sequence_length,
        "max_header_length": max_header_length,
        "sequences": fai_dict,
    }

    with open(output_path, "w") as f:
        json.dump(stats, f, indent=2)


def main():
    args = parse_args()

    fai_dict: dict[str, dict[str, int | str]] = read_fai(args.fai_file)
    fai_dict_assigned: dict[str, dict[str, int | str]] = assign_status(fai_dict, args.max_sequence_length)
    sanitised_dict: dict[str, dict[str, int | str]] = sanitise_headers(fai_dict_assigned, args.max_header_length)

    write_mapping_file(sanitised_dict, f"{args.prefix}_mapping.txt")
    write_new_scaffolds(sanitised_dict, f"{args.prefix}_new_scaffolds.txt")
    write_json_stats(sanitised_dict, args.max_sequence_length, args.max_header_length, f"{args.prefix}_stats.json")


if __name__ == "__main__":
    main()
