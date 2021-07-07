import csv
import os

from argparse import ArgumentParser, FileType


def main():
	parser = ArgumentParser()
	parser.add_argument("source", type=FileType("r"))
	parser.add_argument("target", type=FileType("w"))
	args = parser.parse_args()

	# Get a map of label to ID for all MRO & external terms
	label_id = {}
	with open("index.tsv", "r") as f:
		reader = csv.reader(f, delimiter="\t")
		next(reader)
		next(reader)
		for row in reader:
			label_id[row[1]] = row[0]
	with open("ontology/external.tsv", "r") as f:
		reader = csv.reader(f, delimiter="\t")
		next(reader)
		next(reader)
		for row in reader:
			label_id[row[1]] = row[0]

	# Replace any labels with single quotes with their IDs
	# Only check for cells in "C" template string columns
	rows = []
	reader = csv.DictReader(args.source, delimiter="\t")
	headers = reader.fieldnames
	if None in headers:
		# Extra tab may have been added to the sheet, clean it out
		headers.remove(None)
	robot_strings = next(reader)
	rows.append(robot_strings)
	check_headers = [x for x, y in robot_strings.items() if y.startswith("C ")]
	for row in reader:
		for ch in check_headers:
			value = row[ch]
			if value.startswith("("):
				continue
			if "'" in value:
				term_id = label_id.get(value)
				if not term_id:
					print(f"ERROR: Cannot find ID for '{value}'")
				else:
					row[ch] = term_id
		if None in row:
			# Extra tab may have been added to the sheet, clean it out
			del row[None]
		rows.append(row)

	writer = csv.DictWriter(args.target, delimiter="\t", lineterminator="\n", fieldnames=headers)
	writer.writeheader()
	writer.writerows(rows)


if __name__ == '__main__':
	main()
