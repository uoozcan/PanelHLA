#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Parse SpecHLA's wide genotype table into the pipeline's long format.

SpecHLA writes one row per sample with a column pair per gene:

    # version: IPD-IMGT/HLA 3.38.0
    Sample   HLA_A_1         HLA_A_2         HLA_B_1         ...
    HG00096  A*01:01:01:01   A*29:02:01:01   B*08:01:01:01   ...

Every other caller in this pipeline emits one row per gene, and
modules/aggregation.nf reads only that: it takes field 1 as the gene and fields 2
and 3 as the alleles. Against SpecHLA's layout it found the literal gene "Sample"
on one line and the sample id on the next, matched neither against mv_genes, and
dropped both -- so SpecHLA ran, produced correct calls, and contributed nothing to
the consensus.

The column convention here is the one bin/hla_benchmark.py already uses in
extract_wide_gene_columns(), so the pipeline path and the analysis path read
SpecHLA the same way:

    (?:HLA[_-])?([A-Za-z0-9]+)[_-]([12])$      HLA_A_1 + HLA_A_2 -> gene A

Standard library only: the other bin/parse_*.py scripts import nothing
third-party, and SpecHLA runs natively on Roihu outside any container.

Output:
    # SpecHLA results for {sample}
    Gene    Allele1    Allele2
    A       A*01:01:01:01    A*29:02:01:01
"""

from __future__ import print_function
import argparse
import re
import sys


# Mirrors extract_wide_gene_columns() in bin/hla_benchmark.py. Keep the two in
# step -- tests/test_parse_spechla.py asserts they agree on the same header.
WIDE_COLUMN = re.compile(r"(?:HLA[_-])?([A-Za-z0-9]+)[_-]([12])$")

# What aggregation.nf already treats as a missing allele.
MISSING = {"", "-", "NA", "None", "none", "."}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a SpecHLA hla.result.txt to the pipeline's long format.")
    parser.add_argument("--input", required=True, help="SpecHLA hla.result.txt")
    parser.add_argument("--sample", required=True, help="Sample id for the header")
    parser.add_argument("--output", required=True, help="Long-format TSV to write")
    return parser.parse_args()


def gene_columns(header):
    """Map gene -> (index of allele 1, index of allele 2) from a header row."""
    found = {}
    for index, column in enumerate(header):
        match = WIDE_COLUMN.match(column.strip())
        if not match:
            continue
        gene = match.group(1).strip().upper()
        found.setdefault(gene, {})[match.group(2)] = index
    return dict((gene, (parts["1"], parts["2"]))
                for gene, parts in found.items()
                if "1" in parts and "2" in parts)


def read_rows(path):
    """Header and data rows, ignoring SpecHLA's leading '# version:' comment."""
    header = None
    rows = []
    with open(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
            else:
                rows.append(fields)
    return header, rows


def allele(fields, index):
    value = fields[index].strip() if index < len(fields) else ""
    return "-" if value in MISSING else value


def main():
    args = parse_args()

    header, rows = read_rows(args.input)
    if not header:
        print("SpecHLA file has no header: {}".format(args.input), file=sys.stderr)
        return 1

    columns = gene_columns(header)
    if not columns:
        print("No HLA_<gene>_1/_2 column pairs in {}".format(args.input), file=sys.stderr)
        print("Header was: {}".format("\t".join(header)), file=sys.stderr)
        return 1
    if not rows:
        print("SpecHLA file has a header but no calls: {}".format(args.input), file=sys.stderr)
        return 1

    # SpecHLA writes one row per sample and is run per sample here; take the first.
    fields = rows[0]

    with open(args.output, "w") as out:
        out.write("# SpecHLA results for {}\n".format(args.sample))
        out.write("Gene\tAllele1\tAllele2\n")
        written = 0
        for gene in sorted(columns):
            first, second = columns[gene]
            a1 = allele(fields, first)
            a2 = allele(fields, second)
            if a1 == "-" and a2 == "-":
                continue
            out.write("{}\t{}\t{}\n".format(gene, a1, a2))
            written += 1

    print("SpecHLA parsing complete: {} genes typed for {}".format(written, args.sample))
    return 0


if __name__ == "__main__":
    sys.exit(main())
