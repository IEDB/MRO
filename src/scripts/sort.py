#!/usr/bin/env python3

# This script takes one or more TSV files
# and sorts them by their first column
# using a natural number sort
# and skipping the header rows.

import argparse
import csv
import re


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
# http://stackoverflow.com/a/16090640
    return [
        int(text) if text.isdigit() else text
        for text in re.split(_nsre, s)
    ]

parser = argparse.ArgumentParser(
    description='Sort TSV files')
parser.add_argument('files',
                    type=str,
                    nargs='+',
                    help='TSV files to sort')
args = parser.parse_args()


for filename in args.files:
    header1 = None
    header2 = None
    data = []

    with open(filename, 'r') as f:
        rows = csv.reader(f, delimiter='\t')
        header1 = next(rows)
        header2 = next(rows)
        for row in rows:
            if len(row) == 0:
                continue
            result = [cell.strip() for cell in row]
            if "".join(result) == "":
                continue
            data.append(result)

    data.sort(key=lambda x: natural_sort_key(x[0]))

    with open(filename, 'w') as f:
        writer = csv.writer(f, delimiter='\t',
                            quoting=csv.QUOTE_MINIMAL,
                            lineterminator="\n")
        writer.writerow(header1)
        writer.writerow(header2)
        for row in data:
            while len(row) < len(header1):
                row.append("")
            writer.writerow(row)
