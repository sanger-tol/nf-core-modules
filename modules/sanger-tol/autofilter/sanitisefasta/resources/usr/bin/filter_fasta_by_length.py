#!/usr/bin/env python3

import argparse
import json
import os
import sys
import textwrap
from pathlib import Path
from typing import TypedDict

VERSION = "1.2.0"
DESCRIPTION = f"""
---
Script for filtering a FASTA file by sequence length. By default, sequences shorter than a cutoff value will be removed.
Version = {VERSION}
---

Written by Eerik Aunin (ea10)
Modified by Damon-Lee Pointon (@dp24/@DLBPointon)
Further modified to add JSON logging

"""


class LogStats(TypedDict):
    detailed_changes_truncated: int
    total_sequences: int
    detailed_changes: list[dict]
    filter_mode: str
    sequences_retained: int
    sequences_filtered: int
    has_filtering: bool


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="filter_fasta_by_length",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("in_path", type=str, help="Path to input FASTA file")
    parser.add_argument("--cutoff", type=int, help="Cutoff value for filtering")
    parser.add_argument(
        "-l",
        "--low_pass",
        dest="low_pass",
        action="store_true",
        help="Optional: low pass filtering mode (sequences longer than the cutoff value will be removed)",
    )
    parser.add_argument(
        "--remove_original_fasta",
        action="store_true",
        help="Optional: remove the input FASTA file after creating the filtered FASTA file",
    )
    parser.add_argument(
        "--log_file",
        type=str,
        help="Path to output log file in JSON format",
        default=None,
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION)

    return parser.parse_args(argv)


def print_with_fixed_row_length(seq, max_length):
    """
    Input: 1) a string 2) maximum line length in output
    Output: the input string printed to STDOUT in chunks with fixed maximum line length

    Brought in from Eeriks GPF library

    """
    split_seq = [seq[i : i + max_length] for i in range(0, len(seq), max_length)]
    for line in split_seq:
        print(line)


def file_to_generator(file_path):
    """
    Brought in from Eeriks GPF library
    """
    file_path = Path(file_path)

    if file_path.exists():
        with open(file_path) as f:
            for line in f:
                line = line.rstrip()
                yield line
    else:
        raise FileNotFoundError(f"File not found: {file_path}")


def read_fasta_in_chunks(in_path):
    """
    Input: path to FASTA file
    Output (iterator): string tuples where the first element is a FASTA header and the second element is the corresponding FASTA sequence

    Brought in from Eeriks GPF library

    """
    in_data = file_to_generator(in_path)
    current_seq_header = None
    seq = ""
    for line in in_data:
        if line != "":
            if line[0] == ">":
                if seq != "":
                    yield (current_seq_header, seq)
                seq = ""
                current_seq_header = line[1 : len(line)]
            else:
                seq += line
    if seq != "":
        yield (current_seq_header, seq)


def main(args):
    fasta_path = os.path.abspath(args.in_path)

    # Initialize log statistics
    log_stats: LogStats = {
        "detailed_changes_truncated": 0,
        "total_sequences": 0,
        "detailed_changes": [],
        "filter_mode": "",
        "sequences_retained": 0,
        "sequences_filtered": 0,
        "has_filtering": False,
    }

    if (
        args.cutoff == -1
    ):  # When this script is used as a part of a pipeline, -1 can be assigned as a value for the cutoff to indicate that no filtering should be done
        sys.stderr.write(f"The input FASTA sequences ({fasta_path}) will not be filtered by length\n")
        log_stats["filter_mode"] = "none"

    retained_seq_count = 0
    filtered_seq_count = 0
    total_seq_count = 0

    fasta_data = read_fasta_in_chunks(fasta_path)
    for header, seq in fasta_data:
        total_seq_count += 1

        if args.cutoff == -1:
            print(">" + header)
            print_with_fixed_row_length(seq, 80)
            retained_seq_count += 1
        else:
            seq_len = len(seq)
            seq_len_ok_flag = True

            if args.low_pass:
                if seq_len > args.cutoff:
                    seq_len_ok_flag = False
                    filtered_seq_count += 1
                    log_stats["has_filtering"] = True

                    # Add to detailed changes
                    log_stats["detailed_changes"].append(
                        {
                            "header": header,
                            "change_type": "sequence_filtered",
                            "length": seq_len,
                            "reason": f"Sequence length ({seq_len}) exceeds the length cutoff ({args.cutoff})",
                        }
                    )

                    sys.stderr.write(
                        f"Low pass filtering of FASTA sequences by length: removing sequence {header} from the assembly because its length ({seq_len}) exceeds the length cutoff ({args.cutoff})\n"
                    )
            else:
                if seq_len < args.cutoff:
                    seq_len_ok_flag = False
                    filtered_seq_count += 1
                    log_stats["has_filtering"] = True

                    # Add to detailed changes
                    log_stats["detailed_changes"].append(
                        {
                            "header": header,
                            "change_type": "sequence_filtered",
                            "length": seq_len,
                            "reason": f"Sequence length ({seq_len}) is below the length cutoff ({args.cutoff})",
                        }
                    )

                    sys.stderr.write(
                        f"High pass filtering of FASTA sequences by length: removing sequence {header} from the assembly because its length ({seq_len}) is below the length cutoff ({args.cutoff})\n"
                    )

            if seq_len_ok_flag:
                retained_seq_count += 1
                print(">" + header)
                print_with_fixed_row_length(seq, 80)

    # Update log statistics
    log_stats["total_sequences"] = total_seq_count
    log_stats["sequences_retained"] = retained_seq_count
    log_stats["sequences_filtered"] = filtered_seq_count

    # Limit detailed changes to avoid excessive reporting
    if len(log_stats["detailed_changes"]) > 100:
        log_stats["detailed_changes_truncated"] = True
        log_stats["detailed_changes"] = log_stats["detailed_changes"][:100]
        log_stats["detailed_changes"].append(
            {
                "header": "SUMMARY",
                "change_type": "summary",
                "details": f"Only showing 100 of {filtered_seq_count} filtered sequences.",
            }
        )

    # Write log file if specified
    if args.log_file:
        with open(args.log_file, "w") as f:
            json.dump(log_stats, f, indent=2)

    if args.cutoff != -1 and retained_seq_count == 0:
        sys.stderr.write(
            f"No sequences remain in the FASTA file {fasta_path} after filtering the sequences by length (cutoff: {args.cutoff} bp)\n"
        )
        sys.exit(1)

    if args.remove_original_fasta is True:
        os.remove(fasta_path)


if __name__ == "__main__":
    main(parse_args())
