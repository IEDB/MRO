import csv
import logging

from argparse import ArgumentParser, FileType


def main():
    parser = ArgumentParser()
    parser.add_argument("external", type=FileType("r"))
    parser.add_argument("molecule", type=FileType("r"))
    parser.add_argument("output", type=FileType("w"))
    args = parser.parse_args()

    output = []
    reader = csv.DictReader(args.molecule, delimiter="\t")
    next(reader)
    rows = list(reader)
    label_to_id = {row["Label"]: row["ID"] for row in rows}

    reader = csv.DictReader(args.external, delimiter="\t")
    next(reader)
    for row in reader:
        label_to_id[row["Label"]] = row["ID"]

    for row in rows:
        parent_label = row["Parent"]
        parent_id = label_to_id.get(parent_label)
        if not parent_id:
            logging.error(f"Unable to find ID for parent term '{parent_label}'")
        taxon_label = row["In Taxon"]
        taxon_id = None
        if taxon_label:
            taxon_id = label_to_id.get(taxon_label)
            if not taxon_id:
                logging.error(f"Unable to find ID for taxon term '{taxon_label}'")
        output.append(
            {
                "MRO ID": row["ID"],
                "Label": row["Label"],
                "IEDB Label": row["IEDB Label"],
                "Synonyms": row["Synonyms"],
                "Parent": parent_label,
                "Parent ID": parent_id,
                "In Taxon": taxon_label,
                "In Taxon ID": taxon_id,
            }
        )

    writer = csv.DictWriter(
        args.output,
        fieldnames=[
            "MRO ID",
            "Label",
            "IEDB Label",
            "Synonyms",
            "Parent",
            "Parent ID",
            "In Taxon",
            "In Taxon ID",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(output)


if __name__ == "__main__":
    main()
