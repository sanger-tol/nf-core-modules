#!/usr/bin/env python3

import argparse
import os
import shutil
import sys
from collections.abc import Container, Sized
from pathlib import Path
from typing import Dict, List, Set, TypeAlias

BUSCO_MODES = ("Missing", "Complete", "Duplicated", "Fragmented")

LineageName: TypeAlias = str
BuscoMode: TypeAlias = str
GeneName: TypeAlias = str


def read_full_table(filename: str) -> tuple[LineageName, Dict[BuscoMode, List[GeneName]]]:
    """Read a Busco full_table file and return a tuple of (lineage, genes)
    where `genes` maps modes (cf BUSCO_MODES) to list(gene_ids).
    """
    lineage = None
    genes: Dict[BuscoMode, List[GeneName]] = {s: [] for s in BUSCO_MODES}
    with open(filename) as fh:
        for line in fh:
            if line.startswith("# The lineage dataset is: "):
                lineage = line.split()[5]
            elif line and line[0] != "#":
                t = line[:-1].split("\t")
                gene_id = t[0]
                mode = t[1]
                # Will raise KeyError if the mode isn't recognised
                genes[mode].append(gene_id)
    if lineage is None:
        print(f"Could not parse the lineage name in '{filename}'", file=sys.stderr)
        raise SystemExit(1)
    return (lineage, genes)


def read_full_tables(files: List[str]) -> Dict[LineageName, Dict[BuscoMode, Set[GeneName]]]:
    """Read multiple full_table files and aggregate genes per-lineage.

    Returns `lineage_map` which maps each lineage to a dict of
    mode -> set(gene_ids).
    """
    lineage_map: Dict[LineageName, Dict[BuscoMode, Set[GeneName]]] = {}
    for ft in files:
        lineage, genes = read_full_table(ft)
        if lineage not in lineage_map:
            lineage_map[lineage] = {s: set() for s in BUSCO_MODES}
        for mode, gene_list in genes.items():
            lineage_map[lineage][mode].update(gene_list)
    return lineage_map


def parse_args():
    parser = argparse.ArgumentParser(description="Create a reduced Busco database")
    parser.add_argument("--complete", type=int, help="Number of complete genes to include", default=10)
    parser.add_argument("--duplicated", type=int, help="Number of duplicated genes to include", default=1)
    parser.add_argument("--fragmented", type=int, help="Number of fragmented genes to include", default=2)
    parser.add_argument("--missing", type=int, help="Number of missing genes to include", default=5)
    parser.add_argument("complete_database", help="Path to the complete Busco database")
    parser.add_argument("output_dir", help="Output directory (completely recreated)")
    parser.add_argument(
        "full_table", nargs="+", help="Path(s) to one or more full_table.csv outputs from existing Busco runs"
    )
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")
    return parser.parse_args()


# Collection of helper methods to filter various file formats
# All paths are evaluated relative to `input_dir` and `output_dir`
class BuscoReducer:
    def __init__(self, input_dir: Path, output_dir: Path):
        self.input_dir = input_dir
        self.output_dir = output_dir

    # Filter a TSV file:
    # - `filters` maps column numbers (0-based) to the only values that are allowed in said
    #    column. Rows that have different values in any of those columns are discarded.
    # - `header` is the number of lines to copy as-is.
    def filter_tsv(self, filename: str, filters: Dict[int, Container[str]], header: int = 0):
        input = self.input_dir / filename
        output = self.output_dir / filename
        print("filter_tsv", filename)
        with open(input) as fhi:
            with open(output, "w") as fho:
                for _ in range(header):
                    fho.write(fhi.readline())
                for line in fhi:
                    t = line[:-1].split("\t")
                    if all(t[col] in values for (col, values) in filters.items()):
                        fho.write(line)

    # Filter a FASTA file, keeping only sequences whose identifier starts with one of the
    # values in `values`. The identifier is considered to be the part after the '>' up to
    # the first underscore. That's what matches the gene IDs in the BUSCO full_table files.
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

    # Copy a file as-is
    def copy(self, filename: str):
        input = self.input_dir / filename
        output = self.output_dir / filename
        print("copy", filename)
        with open(input) as fhi:
            with open(output, "w") as fho:
                shutil.copyfileobj(fhi, fho)

    # Write a reduced dataset.cfg file, simply updating the number_of_BUSCOs field
    def write_dataset_cfg(self, filename: str, subset_genes: Sized):
        input = self.input_dir / filename
        output = self.output_dir / filename
        print("write dataset config", filename)
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
    all_full_tables = read_full_tables(args.full_table)

    expected_counts: Dict[BuscoMode, int] = dict(
        zip(
            BUSCO_MODES,
            [
                args.missing,
                args.complete,
                args.duplicated,
                args.fragmented,
            ],
        )
    )

    selected_genes: Dict[LineageName, Dict[GeneName, BuscoMode]] = {}
    all_selected_genes: Set[GeneName] = set()
    all_selected_odb10_genes: Set[GeneName] = set()
    for lineage, gene_map in all_full_tables.items():
        selected_genes[lineage] = {}
        for mode in BUSCO_MODES:
            # Consider genes that have already been selected from other modes
            selected = all_selected_genes.intersection(gene_map[mode])
            if len(selected) < expected_counts[mode]:
                selected.update(sorted(gene_map[mode])[: (expected_counts[mode] - len(selected))])
            for gene in selected:
                selected_genes[lineage][gene] = mode
            all_selected_genes.update(selected)
            if lineage.endswith("_odb10"):
                all_selected_odb10_genes.update(selected)
    all_odb10_lineages = set(filter(lambda lin: lin.endswith("_odb10"), all_full_tables.keys()))

    # Create output directory
    output_dir = Path(args.output_dir)
    input_db = Path(args.complete_database)
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Copy and filter files in the main database directory
    reducer = BuscoReducer(input_db, output_dir)
    reducer.filter_tsv("file_versions.tsv", {0: selected_genes})
    reducer.filter_tsv(
        "info_mappings_all_busco_datasets_odb10.txt", {2: all_selected_odb10_genes, 3: all_odb10_lineages}, 1
    )
    # Then process each lineage
    (output_dir / "lineages").mkdir()
    for lineage, genes in selected_genes.items():
        print("Processing lineage", lineage, genes)
        # Create lineage output directory
        lineage_dir = output_dir / "lineages" / lineage
        lineage_dir.mkdir()
        input_lineage = input_db / "lineages" / lineage
        lineage_reducer = BuscoReducer(input_lineage, lineage_dir)
        # Copy and filter files in the lineage directory
        lineage_reducer.filter_fasta("ancestral", genes)
        lineage_reducer.filter_fasta("ancestral_variants", genes)
        if lineage.endswith("_odb10"):
            lineage_reducer.filter_tsv("lengths_cutoff", {0: genes})
            lineage_reducer.filter_tsv("links_to_ODB10.txt", {0: genes})
        else:
            lineage_reducer.filter_tsv("links_to_ODB12.txt", {0: genes})
        lineage_reducer.filter_fasta("refseq_db.faa", genes)
        lineage_reducer.filter_tsv("scores_cutoff", {0: genes})
        (lineage_dir / "hmms").mkdir()
        for gene in genes:
            lineage_reducer.copy(f"hmms/{gene}.hmm")
        (lineage_dir / "info").mkdir()
        lineage_reducer.filter_tsv("info/ogs.id.info", {1: genes})
        lineage_reducer.copy("info/species.info")
        (lineage_dir / "prfl").mkdir()
        for gene in genes:
            lineage_reducer.copy(f"prfl/{gene}.prfl")

        lineage_reducer.write_dataset_cfg("dataset.cfg", genes)
        write_subset_file(lineage_dir, genes)


if __name__ == "__main__":
    args = parse_args()
    main(args)
