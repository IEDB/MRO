#!/usr/bin/env python3

# Given the chain-sequence.tsv file
# and a copy of ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_prot.fasta
# load the FASTA into a dictionary,
# then for each row of chain-sequence.tsv with an HLA* accession
# replace the sequence with the one from the FASTA file.

import argparse, csv, re

def main():
  parser = argparse.ArgumentParser(
    description='Update HLA sequences')
  parser.add_argument('filename',
      type=str,
      help='TSV file to update')
  parser.add_argument('fasta',
      type=argparse.FileType('r'),
      help='FASTA file to read')
  args = parser.parse_args()

  seqs = {}
  accession = None
  seq = ''
  for line in args.fasta:
    if line.startswith('>'):
      if accession:
        seqs[accession] = seq
      accession = line[5:13]
      seq = ''
    else:
      seq += line.strip()
  seqs[accession] = seq

  rows = []
  with open(args.filename, 'r') as f:
    read_rows = csv.reader(f, delimiter='\t')
    rows = list(read_rows)

  for row in rows:
    if len(row) > 2 and row[3].startswith('HLA'):
      if row[3] in seqs:
        row[4] = seqs[row[3]]
      else:
        print("Missing accession:", row[3])

  with open(args.filename, 'w') as f:
    writer = csv.writer(f, delimiter='\t',
        quoting=csv.QUOTE_MINIMAL,
        lineterminator="\n")
    writer.writerows(rows)

if __name__ == "__main__":
  main()
