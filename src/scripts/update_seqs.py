#!/usr/bin/env python3

# Given the chain-sequence.tsv file
# and one or more FASTA files,
# load the FASTA files into a dictionary,
# then for each row of chain-sequence.tsv with a relevant accession
# replace the sequence with the one from the FASTA file.

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description="Update chain sequences")
    parser.add_argument("filename", type=str, help="TSV file to update")
    parser.add_argument(
        "fasta", type=argparse.FileType("r"), nargs="+", help="FASTA file(s) to read"
    )
    parser.add_argument("-o", "--overwrite", action="store_true", default=False)
    parser.add_argument("-H", "--hla-only", action="store_true", default=False)
    args = parser.parse_args()

    # Read one or more FASTA files into a dictionary
    seqs = {}
    for fasta in args.fasta:
        accession = None
        seq = ""
        for line in fasta:
            if line.startswith(">"):
                if accession:
                    seqs[accession] = seq
                accession = None
                seq = ""

                # Match the accession
                if line.startswith(">HLA:HLA"):
                    accession = line[5:13]
                elif line.startswith(">MHC|SLA"):
                    accession = line[5:13]
                elif line.startswith(">MHC|DLA"):
                    accession = line[5:13]
                elif line.startswith(">IPD-MHC"):
                    accession = line.split(" ")[0][9:]
                else:
                    print("Bad accession:", line)
            else:
                seq += line.strip()
        seqs[accession] = seq

    # Read chain-sequence.tsv
    with open(args.filename, "r") as f:
        read_rows = csv.reader(f, delimiter="\t")
        rows = list(read_rows)

    # Update sequences from FASTAs by matching accession
    for row in rows:
        if len(row) > 2 and row[3] in seqs:
            if args.hla_only and row[2] != "IMGT/HLA":
                continue
            if not row[4] or args.overwrite:
                row[4] = seqs[row[3]]

    # Overwrite chain-sequence.tsv
    with open(args.filename, "w") as f:
        writer = csv.writer(f, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
        writer.writerows(rows)


if __name__ == "__main__":
    main()
