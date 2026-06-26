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


def main():
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

    alns_by_read: dict[str, list] = defaultdict(list)
    with pysam.AlignmentFile(args.input, "rb", threads=args.threads) as bam_in:
        for aln in bam_in:
            if aln.is_unmapped:
                continue
            if args.primary_only and (aln.is_secondary or aln.is_supplementary):
                continue
            if aln.mapping_quality < args.min_mapq:
                continue
            read_id, monomer_idx = parse_read_id(aln.query_name)
            alns_by_read[read_id].append((aln, monomer_idx))

    for read_id in alns_by_read:
        alns = alns_by_read[read_id]
        if alns[0][1] is not None:
            alns.sort(key=lambda x: x[1])
        else:
            alns.sort(key=lambda x: x[0].query_alignment_start)

    with pysam.AlignmentFile(args.input, "rb", threads=args.threads) as bam_in:
        header = bam_in.header.copy()
        with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            with pysam.AlignmentFile(tmp_path, "wb", header=header, threads=args.threads) as bam_out:
                for read_id in sorted(alns_by_read.keys()):
                    alns = alns_by_read[read_id]
                    for frag_idx, (aln, _) in enumerate(alns):
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

            pysam.sort("-o", args.output, tmp_path, "-@", str(args.threads))
            pysam.index(args.output, "-@", str(args.threads))
        finally:
            os.unlink(tmp_path)


if __name__ == "__main__":
    main()
