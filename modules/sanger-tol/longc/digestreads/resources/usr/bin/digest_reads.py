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

"""Digest concatemer reads into monomers at restriction enzyme recognition sites."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from typing import Iterator, TextIO

__version__ = "1.0.0"

# Restriction enzyme recognition sequences (case-insensitive)
ENZYME_SITES = {
    "NlaIII": "CATG",
    "DpnII": "GATC",
    "HindIII": "AAGCTT",
    "MboI": "GATC",
    "Sau3AI": "GATC",
}


def get_site(cutter: str) -> str:
    """Get recognition sequence for enzyme name."""
    if cutter in ENZYME_SITES:
        return ENZYME_SITES[cutter]
    for name, site in ENZYME_SITES.items():
        if name.lower() == cutter.lower():
            return site
    if re.match(r"^[ACGTacgt]+$", cutter):
        return cutter.upper()
    raise ValueError(
        f"Unknown enzyme '{cutter}'. Use one of {list(ENZYME_SITES.keys())} or a recognition sequence (e.g. CATG)"
    )


def split_at_site(seq: str, qual: str, site: str, min_len: int) -> list[tuple[str, str]]:
    """Split sequence and quality at restriction site, return list of (seq, qual) monomers."""
    pattern = re.compile(re.escape(site), re.IGNORECASE)
    monomers = []
    last_end = 0
    for match in pattern.finditer(seq):
        start, end = match.span()
        frag_seq = seq[last_end:start]
        frag_qual = qual[last_end:start] if qual else ""
        if len(frag_seq) >= min_len:
            monomers.append((frag_seq, frag_qual))
        last_end = end
    if last_end < len(seq):
        frag_seq = seq[last_end:]
        frag_qual = qual[last_end:] if qual else ""
        if len(frag_seq) >= min_len:
            monomers.append((frag_seq, frag_qual))
    return monomers


def open_text(path: str) -> TextIO:
    """Open plain, gzipped, or stdin text for reading."""
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


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
    site: str,
    min_len: int,
    out: TextIO,
) -> None:
    """Split one read and write monomer FASTQ records."""
    monomers = split_at_site(seq, qual, site, min_len)
    for idx, (monomer_seq, monomer_qual) in enumerate(monomers):
        monomer_name = f"{name}:{idx}"
        if not monomer_qual:
            monomer_qual = "I" * len(monomer_seq)
        out.write(f"@{monomer_name}\n{monomer_seq}\n+\n{monomer_qual}\n")


def digest_fastq(in_handle: TextIO, out_handle: TextIO, site: str, min_len: int) -> None:
    """Digest all FASTQ records from input handle to output handle."""
    for name, seq, qual in iter_fastq(in_handle):
        emit_monomers(name, seq, qual, site, min_len, out_handle)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Digest concatemers at restriction sites. Reads FASTQ, writes FASTQ monomers."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--cutter", default="NlaIII", help="Restriction enzyme or recognition sequence")
    parser.add_argument("--min-len", type=int, default=1, help="Minimum monomer length (default: 1)")
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        help="Input FASTQ (.fq/.fastq, optionally .gz) or '-' for stdin",
    )
    args = parser.parse_args()

    site = get_site(args.cutter)

    with open_text(args.input) as in_handle:
        digest_fastq(in_handle, sys.stdout, site, args.min_len)


if __name__ == "__main__":
    main()
