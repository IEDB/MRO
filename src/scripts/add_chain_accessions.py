import csv

from argparse import ArgumentParser, FileType


def main():
    parser = ArgumentParser()
    parser.add_argument("iedb", type=FileType("r"))
    parser.add_argument("molecule", type=FileType("r"))
    parser.add_argument("chain", type=FileType("r"))
    parser.add_argument("output", type=FileType("w"))
    args = parser.parse_args()

    chains = {}
    reader = csv.DictReader(args.chain, delimiter="\t")
    next(reader)
    for row in reader:
        if row["Source"] != "IMGT/HLA":
            continue
        chains[row["Label"]] = {"ID": row["ID"], "Accession": row["Accession"]}

    molecules = {}
    reader = csv.DictReader(args.molecule, delimiter="\t")
    next(reader)
    for row in reader:
        chain_i_label = row["Alpha Chain"]
        chain_i_acc = None
        chain_i_id = None
        if chain_i_label:
            chain_i = chains.get(chain_i_label, {})
            if chain_i:
                chain_i_id = chain_i["ID"]
                chain_i_acc = chain_i["Accession"]
        chain_ii_label = row["Beta Chain"]
        chain_ii_acc = None
        chain_ii_id = None
        if chain_ii_label:
            chain_ii = chains.get(chain_i_label, {})
            if chain_ii:
                chain_ii_id = chain_ii["ID"]
                chain_ii_acc = chain_ii["Accession"]
        molecules[row["Label"]] = {
            "chain_i_acc": chain_i_acc,
            "chain_i_id": chain_i_id,
            "chain_ii_acc": chain_ii_acc,
            "chain_ii_id": chain_ii_id
        }

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
