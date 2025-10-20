#!/usr/bin/env python
#
#    Copyright (C) 2025 Genome Research Ltd.
#
#    Author: Jim Downie <jd42@sanger.ac.uk>
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

import gzip

import click
import pyfastx


##
## fastx_format_check borrowed from https://github.com/lmdu/pyfastx/blob/master/pyfastxcli.py
## Author: Lianming Du
## License: MIT
##
def fastx_format_check(fasta):
    if pyfastx.gzip_check(fasta):
        fp = gzip.open(fasta, "rt")
    else:
        fp = open(fasta)

    for line in fp:
        if line.strip():
            break

    fp.close()

    if line[0] == ">":
        return "fasta"
    elif line[0] == "@":
        return "fastq"
    else:
        raise RuntimeError("Error: Input file not FASTA or FASTQ!")


@click.group()
@click.version_option(version="1.0.0", message="%(version)s")
def cli():
    """Main entry point for the tool."""
    pass


@click.command("index")
@click.argument("fasta", type=click.Path(exists=True))
def index(fasta):
    fastx_type = fastx_format_check(fasta)

    if fastx_type == "fasta":
        seq = pyfastx.Fasta(fasta, full_index=False)
    elif fastx_type == "fastq":
        seq = pyfastx.Fastq(fasta, full_index=False)
    else:
        raise RuntimeError("Error: Input file not FASTA or FASTQ!")

    click.echo(len(seq), nl=False)


@click.command("slice")
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("start", type=int)
@click.argument("end", type=int)
def slice(fasta, start, end):
    fastx_type = fastx_format_check(fasta)

    if fastx_type == "fasta":
        seq = pyfastx.Fasta(fasta, full_index=False)
    elif fastx_type == "fastq":
        seq = pyfastx.Fastq(fasta, full_index=False)
    else:
        raise RuntimeError("Error: Input file not FASTA or FASTQ!")

    for i in range(start, end):
        click.echo(seq[i].raw, nl=False)


cli.add_command(index)
cli.add_command(slice)

if __name__ == "__main__":
    cli()
