#!/usr/bin/env python3
#
#    Copyright (C) 2024,2025,2026 Genome Research Ltd.
#
#    Author: Yumi Sims <yy5@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

"""Digest concatemer reads into monomers at Biopython restriction enzyme cut sites.

Cut positions / Bio.Restriction (enzyme.search), not just motif stripping.
"""

from __future__ import annotations

import argparse
import gzip
import math
import sys
from typing import Iterator, TextIO

try:
    from Bio import Restriction
    from Bio.Seq import Seq
except ImportError as exc:
    sys.exit(
        "ERROR: biopython is required for Bio.Restriction digestion "
        f"(import failed: {exc}). Use a container/conda env with biopython."
    )

__version__ = "1.1.0"


def get_enzyme(cutter: str):
    """Resolve a Biopython Restriction enzyme; reject double-cutters (like pore-c-py)."""
    enz = getattr(Restriction, cutter, None)
    if enz is None:
        # Case-insensitive match against Bio.Restriction enzyme names
        for candidate in Restriction.AllEnzymes:
            if candidate.__name__.lower() == cutter.lower():
                enz = candidate
                break
    if enz is None:
        raise ValueError(
            f"Enzyme not found in Bio.Restriction: '{cutter}'. "
            "Use a Biopython enzyme name (e.g. NlaIII, HindIII, DpnII)."
        )
    if enz.cut_twice():
        raise NotImplementedError(f"Enzyme cuts twice, not currently supported: {cutter}")
    return enz


def splits_to_intervals(positions: list[int], length: int) -> list[tuple[int, int]]:
    """Convert 0-based cut positions into monomer [start, end) intervals (pore-c-py)."""
    if len(positions) == 0:
        return [(0, length)]
    prefix = [0] if positions[0] != 0 else []
    suffix = [length] if positions[-1] != length else []
    breaks = prefix + positions + suffix
    return list(zip(breaks[:-1], breaks[1:]))


def cut_intervals(seq: str, enzyme) -> list[tuple[int, int]]:
    """Return monomer intervals using Bio.Restriction cut positions."""
    # Biopython search() is 1-based; pore-c-py converts with x - 1
    cut_points = [x - 1 for x in enzyme.search(Seq(seq))]
    return splits_to_intervals(cut_points, len(seq))


def open_text(path: str) -> TextIO:
    """Open plain, gzipped, or stdin text for reading."""
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def iter_fastq(handle: TextIO) -> Iterator[tuple[str, str, str]]:
    """Yield (name, sequence, quality) records from FASTQ input."""
    while True:
        header = handle.readline()
        if not header:
            return
        if not header.startswith("@"):
            continue
        seq = handle.readline().rstrip("\n")
        plus = handle.readline()
        if not plus:
            return
        qual = handle.readline().rstrip("\n")
        name = header[1:].split()[0]
        yield name, seq, qual


def emit_monomers(
    name: str,
    seq: str,
    qual: str,
    enzyme,
    min_len: int,
    max_monomers: float,
    out: TextIO,
    excluded: TextIO | None,
) -> None:
    """Split one read at enzyme cut sites and write monomer FASTQ records."""
    intervals = cut_intervals(seq, enzyme)
    if len(intervals) > max_monomers:
        if excluded is not None:
            excluded.write(f"{name}\n")
        return

    for idx, (start, end) in enumerate(intervals):
        monomer_seq = seq[start:end]
        if len(monomer_seq) < min_len:
            continue
        monomer_qual = qual[start:end] if qual else "I" * len(monomer_seq)
        if not monomer_qual:
            monomer_qual = "I" * len(monomer_seq)
        out.write(f"@{name}:{idx}\n{monomer_seq}\n+\n{monomer_qual}\n")


def digest_fastq(
    in_handle: TextIO,
    out_handle: TextIO,
    enzyme,
    min_len: int,
    max_monomers: float,
    excluded: TextIO | None,
) -> None:
    """Digest all FASTQ records from input handle to output handle."""
    for name, seq, qual in iter_fastq(in_handle):
        emit_monomers(name, seq, qual, enzyme, min_len, max_monomers, out_handle, excluded)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Digest concatemers at Biopython restriction enzyme cut sites "
            "(pore-c-py style). Reads FASTQ, writes FASTQ monomers."
        )
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--cutter",
        default="NlaIII",
        help="Biopython restriction enzyme name (e.g. NlaIII, HindIII, DpnII)",
    )
    parser.add_argument("--min-len", type=int, default=1, help="Minimum monomer length (default: 1)")
    parser.add_argument(
        "--max-monomers",
        type=int,
        default=None,
        help="Drop reads with more monomers than this (default: no limit; pore-c default 250)",
    )
    parser.add_argument(
        "--excluded-list",
        default=None,
        help="Write names of reads dropped by --max-monomers to this file",
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        help="Input FASTQ (.fq/.fastq, optionally .gz) or '-' for stdin",
    )
    args = parser.parse_args()

    try:
        enzyme = get_enzyme(args.cutter)
    except (ValueError, NotImplementedError) as exc:
        sys.exit(f"ERROR: {exc}")

    max_monomers = float(args.max_monomers) if args.max_monomers is not None else math.inf

    excluded_handle: TextIO | None = None
    try:
        if args.excluded_list:
            excluded_handle = open(args.excluded_list, "w", encoding="utf-8")
        with open_text(args.input) as in_handle:
            digest_fastq(
                in_handle,
                sys.stdout,
                enzyme,
                args.min_len,
                max_monomers,
                excluded_handle,
            )
    except (OSError, ValueError) as exc:
        sys.exit(f"ERROR: {exc}")
    finally:
        if excluded_handle is not None:
            excluded_handle.close()


if __name__ == "__main__":
    main()
