#!/usr/bin/env python3

# This script takes one or more TSV files,
# update the 'synonyms' column
# to include automatically generated synonyms,
# sorts the synonyms,
# and updates the TSV file.

import argparse, csv, re

# http://stackoverflow.com/a/16090640
def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
  return [
    int(text) if text.isdigit() else text.lower()
    for text in re.split(_nsre, s)
  ]

def generate_synonyms(label, synonyms):
  """Given a label and some synonyms, return a set of variations."""
  ss = set()
  for s in synonyms:
    if s.startswith('H2-'):
      ss.add(s.replace('H2-', 'H-2-'))
  for l in [label]:
    ss.add(l.replace('H2-', 'H-2-'))
    ss.add(l.replace('*', ''))
    ss.add(l.replace(':', ''))
    ss.add(l.replace('*','').replace(':', ''))
  return ss

def remove_synonyms(label, synonyms):
  """Given a label and some synonyms,
  return a set of synonyms without the automatic variations."""
  synonyms = set(synonyms)
  synonyms.difference_update(generate_synonyms(label, synonyms))
  synonyms.discard(label)
  synonyms.discard('')
  return synonyms

def update_synonyms(label, synonyms):
  """Given a label and some synonyms, return an updated set of synonyms."""
  synonyms = set(synonyms)
  synonyms.update(generate_synonyms(label, synonyms))
  synonyms.discard(label)
  synonyms.discard('')
  return synonyms

test_cases = [
  ('H2-k class I', set(), {'H-2-k class I'}),
  ('foo', {'H2-Kbm3'}, {'H-2-Kbm3', 'H2-Kbm3'}),
  ('HLA-C*07:01', {'HLA-Cw*0701', 'HLA-Cw*07011'},
    {'HLA-Cw*0701', 'HLA-Cw*07011',
      'HLA-C07:01', 'HLA-C*0701', 'HLA-C0701'}),
  ('HLA-DPA1*01:03/DPB1*03:01', set(),
    {'HLA-DPA101:03/DPB103:01', 'HLA-DPA1*0103/DPB1*0301', 'HLA-DPA10103/DPB10301'})
]

def test_remove_synonyms():
  for case in test_cases:
    (label, synonyms, result) = case
    assert synonyms == remove_synonyms(label, result)

def test_update_synonyms():
  for case in test_cases:
    (label, synonyms, result) = case
    assert result == update_synonyms(label, synonyms)


def main():
  parser = argparse.ArgumentParser(
      description='Add automatically generated synonyms')
  parser.add_argument('infile',
      type=argparse.FileType('r'),
      help='TSV file to process')
  args = parser.parse_args()

  rows = csv.reader(args.infile, delimiter='\t')

  # Find special columns to update
  header = next(rows)
  L = -1
  S = -1
  R = -1
  for i in range(0,len(header)):
    if header[i] in ['IEDB Label', 'displayed_restriction']:
      L = i
    if header[i] in ['Synonyms', 'synonyms']:
      S = i
    if header[i] in ['Restriction Level', 'restriction_level']:
      R = i
  if L == -1:
    raise Exception('Could not find "IEDB Label" column')
  if S == -1:
    raise Exception('Could not find "Synonyms" column')
  if R == -1:
    raise Exception('Could not find "Restriction Level" column')

  print('\t'.join(header))
  for row in rows:
    label = row[L]
    synonyms = row[S].split('|')
    if row[S] != 'A IAO:0000118 SPLIT=|': # ignore special case
      if row[R] in ['complete molecule', 'partial molecule', 'haplotype', 'locus']:
        synonyms = update_synonyms(label, synonyms)
      row[S] = '|'.join(sorted(synonyms, key=natural_sort_key))
    print('\t'.join(row))

if __name__ == "__main__":
  main()

