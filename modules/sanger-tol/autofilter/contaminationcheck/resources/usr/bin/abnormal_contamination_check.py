#!/usr/bin/env python

import argparse
import os.path
import sys
import textwrap
from pathlib import Path

VERSION = "V1.3.0"

DESCRIPTION = """
-------------------------------------
    Abnormal Contamination Check
        Version = {VERSION}
-------------------------------------
Written by James Torrance
Modified by Eerik Aunin
Modified by Damon-Lee Pointon
-------------------------------------

Script for determining if there is
enough contamination found by FCS-GX
to warrant an abnormal contamination
report alarm. Partially based on code
written by James Torrance
-------------------------------------

"""


def parse_args():
    parser = argparse.ArgumentParser(
        prog="Abnormal Contamination Check",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("fai", type=str, help="Path to the assembly FAI file")
    parser.add_argument("summary_path", type=str, help="Path to the tiara summary file")
    parser.add_argument("-q", "--prefix", type=str, help="Output file prefix for the report")
    parser.add_argument("-o", "--output", type=str, help="Path to output file", default="alarm_indicator_file.txt")
    parser.add_argument(
        "-p",
        "--alarm_percentage",
        type=int,
        help="Percentage of putative contaminant sequence in genomic assembly that will trip the alarm",
        default=3,
    )
    parser.add_argument(
        "-l",
        "--alarm_length_removed",
        type=int,
        help="Length of removed sequence is greater than default, greater than this will trip the alarm.",
        default=1e7,
    )
    parser.add_argument(
        "-s", "--alarm_scaff_length", type=int, help="Length of largest scaffold removed to trip alarm.", default=1.8e6
    )
    parser.add_argument(
        "-t",
        "--alarm_scaff_percent_removed",
        type=float,
        help="Percentage of Scaffolds set for removal from assembly to trip the alarm.",
        default=10.0,
    )
    parser.add_argument("-r", "--review_info", type=int, help="Number of REVIEW/INFO to the trigger alarm", default=0)
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    return parser.parse_args()


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


def get_sequence_lengths(assembly_fai):
    """
    Read FAI file and return a dictionary of sequence lengths
    """
    seq_lengths_dict = dict()
    with open(assembly_fai) as fai_file:
        for line in fai_file:
            split_line = line.split("\t")
            seq_name = split_line[0]
            seq_len = int(split_line[1])
            seq_lengths_dict[seq_name] = dict()
            seq_lengths_dict[seq_name]["seq_len"] = seq_len
    return seq_lengths_dict


def load_fcs_gx_results(seq_dict, fcs_gx_and_tiara_summary_path):
    """
    Loads FCS-GX actions from the FCS-GX and Tiara results summary file, adds them to the dictionary that contains sequence lengths
    """
    fcs_gx_and_tiara_summary_data = file_to_generator(fcs_gx_and_tiara_summary_path)
    # Skip the header line containing column names
    next(fcs_gx_and_tiara_summary_data)

    for line in fcs_gx_and_tiara_summary_data:
        split_line = line.split(",")
        assert len(split_line) == 5
        seq_name = split_line[0]
        fcs_gx_action = split_line[1]
        seq_dict[seq_name]["fcs_gx_action"] = fcs_gx_action

    return seq_dict


def main():
    args = parse_args()
    if not os.path.isfile(args.summary_path):
        sys.stderr.write(
            f"The FCS-GX and Tiara results file was not found at the expected location ({args.summary_path})\n"
        )
        sys.exit(1)

    if not os.path.isfile(args.fai):
        sys.stderr.write(f"The assembly FASTA file was not found at the expected location ({args.fai})\n")
        sys.exit(1)

    seq_dict = get_sequence_lengths(args.fai)
    seq_dict = load_fcs_gx_results(seq_dict, args.summary_path)

    total_assembly_length = 0
    lengths_removed = list()
    scaffolds_removed = 0
    scaffold_count = len(seq_dict)
    review_info = 0

    for seq_name in seq_dict:
        seq_len = seq_dict[seq_name]["seq_len"]
        if seq_dict[seq_name]["fcs_gx_action"] == "EXCLUDE":
            lengths_removed.append(seq_len)
            scaffolds_removed += 1
        if seq_dict[seq_name]["fcs_gx_action"] in ["REVIEW", "INFO"]:
            review_info += 1
        total_assembly_length += seq_len

    alarm_threshold_for_parameter = {
        "TOTAL_LENGTH_REMOVED": args.alarm_length_removed,
        "PERCENTAGE_LENGTH_REMOVED": args.alarm_percentage,
        "LARGEST_SCAFFOLD_REMOVED": args.alarm_scaff_length,
        "PERCENTAGE_SCAFFOLDS_REMOVED": args.alarm_scaff_percent_removed,
        "REVIEW_OR_INFO": args.review_info,
    }

    report_dict = {
        "TOTAL_LENGTH_REMOVED": sum(lengths_removed),
        "PERCENTAGE_LENGTH_REMOVED": 100 * sum(lengths_removed) / total_assembly_length,
        "LARGEST_SCAFFOLD_REMOVED": max(lengths_removed, default=0),
        "SCAFFOLDS_REMOVED": scaffolds_removed,
        "PERCENTAGE_SCAFFOLDS_REMOVED": 100 * scaffolds_removed / scaffold_count,
        "REVIEW_OR_INFO": review_info,
    }

    with open(f"{args.prefix}_raw_report.txt", "a") as f:
        for x, y in report_dict.items():
            print(f"{x}: {y}", file=f)

    alarm_list = []
    for param in alarm_threshold_for_parameter:
        param_value = report_dict[param]
        alarm_threshold = alarm_threshold_for_parameter[param]

        # IF CONTAMINATING SEQ FOUND FILL FILE WITH ABNORMAL CONTAM
        if param_value > alarm_threshold_for_parameter[param]:
            alarm_list.append(
                f"YES_ABNORMAL_CONTAMINATION: Stage 1 decon for {args.fai}: {param} == {param_value} : Alarm threshold == {alarm_threshold}\n"
            )
        else:
            alarm_list.append(
                f"NO_ABNORMAL_CONTAMINATION: Stage 1 decon for {args.fai}: {param} == {param_value} : Alarm threshold == {alarm_threshold}\n"
            )

    # Separated out to ensure that the file is written in one go and doesn't confuse Nextflow
    with open(args.output, "w") as f:
        f.write("".join(alarm_list))


if __name__ == "__main__":
    main()
