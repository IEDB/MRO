import csv
import logging

from argparse import ArgumentParser, FileType
from collections import defaultdict


def main():
	parser = ArgumentParser()
	parser.add_argument("index", type=FileType("r"))
	parser.add_argument("external", type=FileType("r"))
	parser.add_argument("molecule", type=FileType("r"))
	parser.add_argument("output", type=FileType("w"))
	args = parser.parse_args()

	label_to_id = {}
	reader = csv.DictReader(args.index, delimiter="\t")
	next(reader)
	for row in reader:
		label_to_id[row["Label"]] = row["ID"]

	reader = csv.DictReader(args.external, delimiter="\t")
	next(reader)
	for row in reader:
		label_to_id[row["Label"]] = row["ID"]

	output = []
	reader = csv.DictReader(args.molecule, delimiter="\t")
	next(reader)
	for row in reader:
		term_label = row["Label"]
		term_id = label_to_id.get(term_label)
		if not term_id:
			logging.error(f"Unable to find ID for '{term_label}'")
			continue
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
		output.append({"MRO ID": term_id, "Label": term_label, "IEDB Label": row["IEDB Label"], "Synonyms": row["Synonyms"], "Parent": parent_id, "In Taxon": taxon_label, "In Taxon ID": taxon_id})

	writer = csv.DictWriter(args.output, fieldnames=["MRO ID", "Label", "IEDB Label", "Synonyms", "Parent", "In Taxon", "In Taxon ID"], delimiter="\t", lineterminator="\n")
	writer.writeheader()
	writer.writerows(output)


if __name__ == '__main__':
	main()

