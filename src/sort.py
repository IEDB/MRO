#!/usr/bin/env python3

# This script takes one or more TSV files
# and sorts them by their first column
# using a natural number sort
# and skipping the header rows.

import argparse, csv, re

# http://stackoverflow.com/a/16090640
def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
  return [
    int(text) if text.isdigit() else text.lower()
    for text in re.split(_nsre, s)
  ]

# Parse arguments

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
      data.append(row)

  data.sort(key=lambda x: natural_sort_key(x[0]))

  with open(filename, 'w') as f:
    writer = csv.writer(f, delimiter='\t', lineterminator="\n")
    writer.writerow(header1)
    writer.writerow(header2)
    writer.writerows(data)

