#!/usr/bin/env python3

# Given the chain.tsv file
# and one or more FASTA files,
# load the FASTA files into a dictionary,
# then for each row of chain.tsv with a relevant accession
# replace the sequence with the one from the FASTA file.

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description="Update chain sequences")
    parser.add_argument("chain", type=str, help="chain TSV file to update")
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

    # Read chain.tsv
    with open(args.chain, "r") as f:
        read_rows = csv.DictReader(f, delimiter="\t")
        header = read_rows.fieldnames
        robot_header = next(read_rows)
        rows = list(read_rows)

    # Update sequences from FASTAs by matching accession
    new_rows = [robot_header]
    for row in rows:
        acc = row["Accession"]
        if acc and acc in seqs:
            if args.hla_only and row["Source"] != "IMGT/HLA":
                continue
            if not row["Sequence"] or args.overwrite:
                row["Sequence"] = seqs[acc]
        new_rows.append(row)

    # Overwrite chain.tsv
    with open(args.chain, "w") as f:
        writer = csv.DictWriter(
            f, fieldnames=header, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(new_rows)


if __name__ == "__main__":
    main()
