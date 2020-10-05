#!/usr/bin/env python3

# This script takes all TSV tables and looks for new labels that do not exist in index.tsv
# These labels are added to index.tsv and assigned new IDs
# It also looks at iedb.tsv and assigns new numeric IDs to newly-added terms.

import csv

from argparse import ArgumentParser


def main():
	p = ArgumentParser()
	p.add_argument("index")
	p.add_argument("iedb")
	p.add_argument("templates_dir")
	args = p.parse_args()

	# Update iedb.tsv by assigning numeric IDs to any new terms
	# These should always exist in other sheets
	
	iedb_id = 0
	with open(args.iedb, "r") as f:
		iedb_rows = []
		reader = csv.DictReader(f, delimiter="\t")
		iedb_headers = reader.fieldnames
		# Skip ROBOT template
		iedb_rows.append(next(reader))
		for row in reader:
			iedb_rows.append(row)
			num_id = row["IEDB ID"]
			if num_id and num_id.strip() != "":
				num_id = int(num_id)
				if num_id > iedb_id:
					iedb_id = num_id
	added_iedb = 0
	with open(args.iedb, "w") as f:
		writer = csv.DictWriter(f, delimiter="\t", lineterminator="\n", fieldnames=iedb_headers)
		writer.writeheader()
		for row in iedb_rows:
			num_id = row["IEDB ID"]
			if not num_id or num_id.strip() == "":
				added_iedb += 1
				iedb_id += 1
				row["IEDB ID"] = iedb_id
			writer.writerow(row)

	if added_iedb:
		print(f"\nAdded {added_iedb} new IDs to {args.iedb}\n")

	# Read in index.tsv to get current highest-ID and current rows

	existing_labels = []
	rows = []
	mro_id = 0
	with open(args.index, "r") as f:
		reader = csv.DictReader(f, delimiter="\t")
		headers = reader.fieldnames
		# Skip ROBOT template strings
		rows.append(next(reader))
		for row in reader:
			rows.append(row)
			if row["Type"] != "owl:Class":
				continue
			existing_labels.append(row["Label"])
			num_id = int(row["ID"].split(":")[1].lstrip("0"))
			if num_id > mro_id:
				mro_id = num_id

	# Read in tables in template directory to find labels not in index

	new_labels = []
	for tbl in ["chain", "genetic-locus", "haplotype", "haplotype-molecule", "molecule", "mutant-molecule", "serotype", "serotype-molecule"]:
		with open(f"{args.templates_dir}/{tbl}.tsv", "r") as f:
			reader = csv.DictReader(f, delimiter="\t")
			# Skip ROBOT template strings
			next(reader)
			for row in reader:
				label = row["Label"]
				if label not in existing_labels:
					new_labels.append(label)

	if not new_labels:
		print("\nNo new labels to add!\n")
		return

	# Add all new terms to index.tsv
	
	print(f"\nAdding {len(new_labels)} new term(s) to index\n")

	for nl in new_labels:
		mro_id += 1
		rows.append({"ID": f"MRO:{mro_id:07}", "Label": nl, "Type": "owl:Class"})

	with open(args.index, "w") as f:
		writer = csv.DictWriter(f, delimiter="\t", lineterminator="\n", fieldnames=headers)
		writer.writeheader()
		writer.writerows(rows)


if __name__ == '__main__':
	main()
