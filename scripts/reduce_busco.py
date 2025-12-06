#!/usr/bin/env python3

import argparse
import os
import shutil
from collections.abc import Container, Sized
from pathlib import Path
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


class BuscoReducer:
    def __init__(self, input_dir: Path, output_dir: Path):
        self.input_dir = input_dir
        self.output_dir = output_dir

    def filter_tsv(self, filename: str, filters: Dict[int, Container[str]], header: int):
        input = self.input_dir / filename
        output = self.output_dir / filename
        print("filter_tsv", filename)
        with open(input) as fhi:
            with open(output, "w") as fho:
                for i, line in enumerate(fhi):
                    if i < header:
                        fho.write(line)
                    else:
                        t = line[:-1].split("\t")
                        if all(t[col] in values for (col, values) in filters.items()):
                            fho.write(line)

    def filter_fasta(self, filename: str, values: Container[str]):
        input = self.input_dir / filename
        output = self.output_dir / filename
        print("filter_fasta", filename)
        copy_block = False
        with open(input) as fhi:
            with open(output, "w") as fho:
                for line in fhi:
                    if line.startswith(">"):
                        copy_block = line[1:-1].split("_")[0] in values
                    if copy_block:
                        fho.write(line)

    def copy(self, filename: str):
        input = self.input_dir / filename
        output = self.output_dir / filename
        print("copy", filename)
        with open(input) as fhi:
            with open(output, "w") as fho:
                shutil.copyfileobj(fhi, fho)

    def write_dataset_cfg(self, filename: str, subset_genes: Sized):
        input = self.input_dir / filename
        output = self.output_dir / filename
        with open(input) as fhi:
            with open(output, "w") as fho:
                for line in fhi:
                    if line.startswith("number_of_BUSCOs="):
                        print(f"number_of_BUSCOs={len(subset_genes)}", file=fho)
                    else:
                        fho.write(line)


def write_subset_file(output_dir: Path, subset_genes: Dict[str, str]):
    with open(output_dir / "SUBSET", "w") as fh:
        for gene, mode in subset_genes.items():
            print(gene, mode, file=fh)


def main(args):
    print(args)
    (lineage, genes) = read_full_table(args.full_table)

    subset_genes: Dict[str, str] = {}
    for mode, count in [
        ("Complete", args.complete),
        ("Duplicated", args.duplicated),
        ("Fragmented", args.fragmented),
        ("Missing", args.missing),
    ]:
        if count > len(genes[mode]):
            print(f"Requested {count} {mode} genes but only {len(genes[mode])} available")
        for gene in sorted(genes[mode])[:count]:
            subset_genes[gene] = mode
    print(f"Selected {len(subset_genes)} genes for reduced database")
    print(subset_genes)

    output_dir = Path(args.output_dir)
    input_db = Path(args.complete_database)
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    reducer = BuscoReducer(input_db, output_dir)
    reducer.filter_tsv("file_versions.tsv", {0: [lineage]}, 0)
    reducer.filter_tsv("info_mappings_all_busco_datasets_odb10.txt", {2: subset_genes, 3: [lineage]}, 0)
    (output_dir / "lineages").mkdir()
    lineage_dir = output_dir / "lineages" / lineage
    lineage_dir.mkdir()
    input_lineage = input_db / "lineages" / lineage
    lineage_reducer = BuscoReducer(input_lineage, lineage_dir)
    lineage_reducer.filter_fasta("ancestral", subset_genes)
    lineage_reducer.filter_fasta("ancestral_variants", subset_genes)
    if lineage.endswith("_odb10"):
        lineage_reducer.filter_tsv("lengths_cutoff", {0: subset_genes}, 0)
        lineage_reducer.filter_tsv("links_to_ODB10.txt", {0: subset_genes}, 0)
    else:
        lineage_reducer.filter_tsv("links_to_ODB12.txt", {0: subset_genes}, 0)
    lineage_reducer.filter_fasta("refseq_db.faa", subset_genes)
    lineage_reducer.filter_tsv("scores_cutoff", {0: subset_genes}, 0)
    (lineage_dir / "hmms").mkdir()
    for gene in subset_genes:
        lineage_reducer.copy(f"hmms/{gene}.hmm")
    (lineage_dir / "info").mkdir()
    lineage_reducer.filter_tsv("info/ogs.id.info", {1: subset_genes}, 0)
    lineage_reducer.copy("info/species.info")
    (lineage_dir / "prfl").mkdir()
    for gene in subset_genes:
        lineage_reducer.copy(f"prfl/{gene}.prfl")

    lineage_reducer.write_dataset_cfg("dataset.cfg", subset_genes)

    write_subset_file(output_dir, subset_genes)


if __name__ == "__main__":
    args = parse_args()
    main(args)
