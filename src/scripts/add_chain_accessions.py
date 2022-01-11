import csv

from argparse import ArgumentParser, FileType


def main():
    parser = ArgumentParser()
    parser.add_argument("index", type=FileType("r"))
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

    label_to_id = {}
    reader = csv.DictReader(args.index, delimiter="\t")
    next(reader)
    for row in reader:
        if row["Label"] not in chain_accessions:
            continue
        label_to_id[row["Label"]] = row["ID"]

    molecules = {}
    reader = csv.DictReader(args.molecule, delimiter="\t")
    next(reader)
    for row in reader:
        chain_i_label = row["Alpha Chain"]
        chain_i_acc = None
        chain_i_id = None
        if chain_i_label:
            chain_i_acc = chain_accessions.get(chain_i_label)
            chain_i_id = label_to_id.get(chain_i_label)
        chain_ii_label = row["Beta Chain"]
        chain_ii_acc = None
        chain_ii_id = None
        if chain_ii_label:
            chain_ii_acc = chain_accessions.get(chain_ii_label)
            chain_ii_id = label_to_id.get(chain_ii_label)
        molecules[row["Label"]] = {"chain_i_acc": chain_i_acc, "chain_i_id": chain_i_id, "chain_ii_acc": chain_ii_acc, "chain_ii_id": chain_ii_id}

    rows = []
    reader = csv.DictReader(args.iedb, delimiter="\t")
    row = next(reader)
    row["Chain I Accession"] = "A MRO:has-chain-i-accession"
    row["Chain I MRO ID"] = "A MRO:has-chain-i-id"
    row["Chain II Accession"] = "A MRO:has-chain-ii-accession"
    row["Chain II MRO ID"] = "A MRO:has-chain-ii-id"
    rows.append(row)
    for row in reader:
        label = row["Label"]
        if label in molecules:
            row["Chain I Accession"] = molecules[label]["chain_i_acc"]
            row["Chain I MRO ID"] = molecules[label]["chain_i_id"]
            row["Chain II Accession"] = molecules[label]["chain_ii_acc"]
            row["Chain II MRO ID"] = molecules[label]["chain_ii_id"]
        rows.append(row)

    writer = csv.DictWriter(
        args.output,
        fieldnames=[
            "Label",
            "IEDB ID",
            "Locus",
            "Chain I Source ID",
            "Chain I Accession",
            "Chain I MRO ID",
            "Chain II Source ID",
            "Chain II Accession",
            "Chain II MRO ID",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
