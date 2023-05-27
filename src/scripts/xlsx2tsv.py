#!/usr/bin/env python3

import argparse, os, csv
from openpyxl import load_workbook

parser = argparse.ArgumentParser(
    description='Read TSV from Excel')
parser.add_argument('input',
    type=str,
    help='Excel file to read')
parser.add_argument('sheet',
    type=str,
    help='Sheet name to read')
args = parser.parse_args()

wb = load_workbook(args.input, read_only=True)
ws = wb[args.sheet]

for row in ws:
  values = []
  for cell in row:
    if cell.value is None:
      values.append('')
    else:
      values.append(str(cell.value))
  print('\t'.join(values))

