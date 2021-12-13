import csv

from argparse import ArgumentParser, FileType


def main():
    parser = ArgumentParser()
    parser.add_argument("iedb", type=FileType("r"))
    parser.add_argument("molecule", type=FileType("r"))
    parser.add_argument("chain_sequence", type=FileType("r"))
    parser.add_argument("output", type=FileType("w"))
    args = parser.parse_args()

    chain_accessions = {}
    reader = csv.DictReader(args.chain_sequence, delimiter="\t")
    next(reader)
    for row in reader:
        if row["Source"] != "IMGT/HLA":
            continue
        chain_accessions[row["Label"]] = row["Accession"]

    molecules = {}
    reader = csv.DictReader(args.molecule, delimiter="\t")
    next(reader)
    for row in reader:
        chain_i_label = row["Alpha Chain"]
        chain_i_acc = None
        if chain_i_label:
            chain_i_acc = chain_accessions.get(chain_i_label)
        chain_ii_label = row["Beta Chain"]
        chain_ii_acc = None
        if chain_ii_label:
            chain_ii_acc = chain_accessions.get(chain_ii_label)
        molecules[row["Label"]] = {"chain_i": chain_i_acc, "chain_ii": chain_ii_acc}

    rows = []
    reader = csv.DictReader(args.iedb, delimiter="\t")
    row = next(reader)
    row["Chain I Accession"] = "A MRO:has-chain-i-accession"
    row["Chain II Accession"] = "A MRO:has-chain-ii-accession"
    rows.append(row)
    for row in reader:
        label = row["Label"]
        if label in molecules:
            row["Chain I Accession"] = molecules[label]["chain_i"]
            row["Chain II Accession"] = molecules[label]["chain_ii"]
        rows.append(row)

    writer = csv.DictWriter(args.output, fieldnames=["Label", "IEDB ID", "Locus", "Chain I Source ID", "Chain I Accession", "Chain II Source ID", "Chain II Accession"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == '__main__':
    main()
