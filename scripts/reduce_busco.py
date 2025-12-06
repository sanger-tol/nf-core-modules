#!/usr/bin/env python3

import argparse
from collections.abc import Container, Sized
import os
from pathlib import Path
import shutil
from typing import Dict, List


def read_full_table(filename):
    """Read a Busco full_table file and return a tuple of (lineage, genes)
    where `genes` is a maps modes -> list(gene_ids).
    """
    lineage = None
    genes: Dict[str, List[str]] = {s: [] for s in ("Missing", "Complete", "Duplicated", "Fragmented")}
    with open(filename) as fh:
        for line in fh:
            if line.startswith("# The lineage dataset is: "):
                lineage = line.split()[5]
            elif line and line[0] != "#":
                t = line[:-1].split("\t")
                genes[t[1]].append(t[0])
    assert lineage
    return (lineage, genes)


def parse_args():
    parser = argparse.ArgumentParser(description="Create a reduced Busco database")
    parser.add_argument("--complete", type=int, help="Number of complete genes to include", default=10)
    parser.add_argument("--duplicated", type=int, help="Number of duplicated genes to include", default=1)
    parser.add_argument("--fragmented", type=int, help="Number of fragmented genes to include", default=2)
    parser.add_argument("--missing", type=int, help="Number of missing genes to include", default=5)
    parser.add_argument("complete_database", help="Path to the complete Busco database")
    parser.add_argument("full_table", help="Path to the full_table.csv output from an existing Busco run")
    parser.add_argument("output_dir", help="Output directory (completely recreated)")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")
    return parser.parse_args()

def filter_tsv(filename: Path, output: Path, filters: Dict[int, Container[str]], header: int):
    print("filter_tsv", filename.name)
    with open(filename) as fhi:
        with open(output, "w") as fho:
            for (i, line) in enumerate(fhi):
                if i < header:
                    fho.write(line)
                else:
                    t = line[:-1].split("\t")
                    if all(t[col] in values for (col, values) in filters.items()):
                        fho.write(line)

def filter_fasta(filename: Path, output: Path, values: Container[str]):
    print("filter_fasta", filename.name)
    copy = False
    with open(filename) as fhi:
        with open(output, "w") as fho:
            for line in fhi:
                if line.startswith(">"):
                    copy = line[1:-1].split("_") in values
                if copy:
                    fho.write(line)

def copy(filename, output):
    print("copy", filename.name)
    with open(filename) as fhi:
        with open(output, "w") as fho:
            for line in fhi:
                fho.write(line)

def main(args):
    print(args)
    (lineage, genes) = read_full_table(args.full_table)
    subset_genes = sorted(genes["Complete"])[:args.complete] \
        + sorted(genes["Duplicated"])[:args.duplicated] \
        + sorted(genes["Fragmented"])[:args.fragmented] \
        + sorted(genes["Missing"])[:args.missing]
    output_dir = Path(args.output_dir)
    input_db = Path(args.complete_database)
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    output_dir.mkdir()
    filter_tsv(input_db / "file_versions.tsv", output_dir / "file_versions.tsv", {0: [lineage]}, 0)
    filter_tsv(input_db / "info_mappings_all_busco_datasets_odb10.txt", output_dir / "info_mappings_all_busco_datasets_odb10.txt", {2: subset_genes, 3: [lineage]}, 0)
    (output_dir / "lineages").mkdir()
    lineage_dir = output_dir / "lineages" / lineage
    lineage_dir.mkdir()
    input_lineage = input_db / "lineages" / lineage
    filter_fasta(input_lineage / "ancestral", lineage_dir / "ancestral", subset_genes)
    filter_fasta(input_lineage / "ancestral_variants", lineage_dir / "ancestral_variants", subset_genes)
    if lineage.endswith("_odb10"):
        filter_tsv(input_lineage / "lengths_cutoff", lineage_dir / "lengths_cutoff", {0: subset_genes}, 0)
        filter_tsv(input_lineage / "links_to_ODB10.txt", lineage_dir / "links_to_ODB10.txt", {0: subset_genes}, 0)
    else:
        filter_tsv(input_lineage / "links_to_ODB12.txt", lineage_dir / "links_to_ODB12.txt", {0: subset_genes}, 0)
    filter_fasta(input_lineage / "refseq_db.faa", lineage_dir / "refseq_db.faa", subset_genes)
    filter_tsv(input_lineage / "scores_cutoff", lineage_dir / "scores_cutoff", {0: subset_genes}, 0)
    (lineage_dir / "hmms").mkdir()
    for gene in subset_genes:
        copy(input_lineage / "hmms" / f"{gene}.hmm", lineage_dir / "hmms" / f"{gene}.hmm")
    (lineage_dir / "info").mkdir()
    filter_tsv(input_lineage / "info" / "ogs.id.info", lineage_dir / "info" / "ogs.id.info", {1: subset_genes}, 0)
    copy(input_lineage / "info" / "species.info", lineage_dir / "info" / "species.info")
    (lineage_dir / "prfl").mkdir()
    for gene in subset_genes:
        copy(input_lineage / "prfl" / f"{gene}.prfl", lineage_dir / "prfl" / f"{gene}.prfl")

    with open(output_dir / "SUBSET", "w") as fh:
        for gene in subset_genes:
            mode = [k for (k,l) in genes.items() if gene in l][0]
            print(gene, mode, file=fh)

    with open(input_lineage / "dataset.cfg") as fhi:
        with open(lineage_dir / "dataset.cfg", "w") as fho:
            for line in fhi:
                if line.startswith("number_of_BUSCOs="):
                    print(f"number_of_BUSCOs={len(subset_genes)}", file=fho)
                else:
                    fho.write(line)


if __name__ == "__main__":
    args = parse_args()
    main(args)
