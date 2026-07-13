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

"""Annotate fragments in aligned BAM: group by read, filter, order, assign fragment IDs."""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from collections import defaultdict

try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam required. Install with: pip install pysam")

__version__ = "1.0.0"


def parse_read_id(qname: str) -> tuple[str, int | None]:
    """Extract read_id and monomer_idx from query name."""
    if ":" in qname:
        parts = qname.rsplit(":", 1)
        if len(parts) == 2 and parts[1].isdigit():
            return parts[0], int(parts[1])
    return qname, None


def sort_read_alignments(alns: list[tuple[pysam.AlignedSegment, int | None]]) -> None:
    """Order monomers by index in read name, else by alignment start."""
    if alns[0][1] is not None:
        alns.sort(key=lambda x: x[1])
    else:
        alns.sort(key=lambda x: x[0].query_alignment_start)


def collect_alignments(
    bam_path: str,
    min_mapq: int,
    primary_only: bool,
    threads: int,
) -> dict[str, list[tuple[pysam.AlignedSegment, int | None]]]:
    """Read BAM, group alignments by concatemer read, and order monomers."""
    alns_by_read: dict[str, list[tuple[pysam.AlignedSegment, int | None]]] = defaultdict(list)
    with pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam_in:
        for aln in bam_in:
            if aln.is_unmapped:
                continue
            if primary_only and (aln.is_secondary or aln.is_supplementary):
                continue
            if aln.mapping_quality < min_mapq:
                continue
            read_id, monomer_idx = parse_read_id(aln.query_name)
            alns_by_read[read_id].append((aln, monomer_idx))

    for alns in alns_by_read.values():
        sort_read_alignments(alns)
    return alns_by_read


def write_annotated_bam(
    bam_path: str,
    output_path: str,
    alns_by_read: dict[str, list[tuple[pysam.AlignedSegment, int | None]]],
    threads: int,
) -> None:
    """Write tagged alignments to a temp BAM, then sort and index the output."""
    with pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam_in:
        header = bam_in.header.copy()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = os.path.join(tmpdir, "unsorted.bam")
            with pysam.AlignmentFile(tmp_path, "wb", header=header, threads=threads) as bam_out:
                for read_id in sorted(alns_by_read):
                    for frag_idx, (aln, _) in enumerate(alns_by_read[read_id]):
                        ref_id = aln.reference_id
                        ref_name = bam_in.get_reference_name(ref_id) if ref_id >= 0 else "."
                        ref_start = aln.reference_start + 1
                        ref_end = aln.reference_end
                        frag_id = f"{ref_name}:{ref_start}-{ref_end}"

                        new_aln = pysam.AlignedSegment.from_dict(aln.to_dict(), header)
                        new_aln.set_tag("FI", frag_idx, "i")
                        new_aln.set_tag("FD", frag_id, "Z")
                        new_aln.set_tag("BX", read_id, "Z")
                        bam_out.write(new_aln)

            pysam.sort("-o", output_path, tmp_path, "-@", str(threads))
            pysam.index(output_path, "-@", str(threads))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate fragments in BAM: group by read, filter, order, assign fragment IDs"
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--input", "-i", required=True, help="Input BAM")
    parser.add_argument("--output", "-o", required=True, help="Output annotated BAM")
    parser.add_argument("--min-mapq", type=int, default=0, help="Minimum mapping quality [0]")
    parser.add_argument("--primary-only", action="store_true", default=True, help="Keep only primary alignments")
    parser.add_argument("--threads", type=int, default=1, help="Threads for BAM I/O")
    args = parser.parse_args()

    try:
        alns_by_read = collect_alignments(args.input, args.min_mapq, args.primary_only, args.threads)
        write_annotated_bam(args.input, args.output, alns_by_read, args.threads)
    except (OSError, ValueError) as exc:
        sys.exit(f"ERROR: {exc}")
    except pysam.utils.SamtoolsError as exc:
        sys.exit(f"ERROR: samtools failed: {exc}")


if __name__ == "__main__":
    main()
