#!/usr/bin/env python3

# This script takes the CSV file with results
# from the mhc_allele_restriction.rq SPARQL query
# and reformats the rows to align with
# the mhc_allele_restriction table in the current IEDB.
# Some of these manipulations could be done in SPARQL,
# but that query is already ugly enough.

import argparse, csv, re
from synonyms import remove_synonyms
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

# http://stackoverflow.com/a/16090640
def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
  return [
    int(text) if text.isdigit() else text.lower()
    for text in re.split(_nsre, s)
  ]

# Parse arguments

parser = argparse.ArgumentParser(
    description='Clean SPARQL results to align with mhc_allele_restriction table')
parser.add_argument('alleles',
    type=argparse.FileType('r'),
    help='read query result CSV')
args = parser.parse_args()

organisms = {
  'bonobo': 'bonobo (Pan paniscus)',
  'cattle': 'cattle (Bos taurus)',
  'chicken': 'chicken (Gallus gallus)',
  'chimpanzee': 'chimpanzee (Pan troglodytes)',
  'cotton-top tamarin': 'cotton-top tamarin (Saguinus oedipus)',
  'crab-eating macaque': 'crab-eating macaque (Macaca fascicularis)',
  'dog': 'dog (Canis lupus familiaris)',
  'duck': 'duck (Anas platyrhynchos)',
  'gorilla': 'gorilla (Gorilla gorilla)',
  'horse': 'horse (Equus caballus)',
  'human': 'human (Homo sapiens)',
  'marmoset': 'marmoset (Callithrix jacchus)',
  'mouse': 'mouse (Mus musculus)',
  'organism': 'organism (all species)',
  'pig': 'pig (Sus scrofa)',
  'rat': 'rat (Rattus norvegicus)',
  'rhesus macaque': 'rhesus macaque (Macaca mulatta)',
  'sheep': 'sheep (Ovis aries)'
}
codes = {
  '1': 'organism',
  '8839': 'duck',
  '9031': 'chicken',
  '9483': 'marmoset',
  '9490': 'Saoe',
  '9544': 'Mamu',
  '9541': 'Mafa',
  '9593': 'Gogo',
  '9597': 'Papa',
  '9598': 'Patr',
  '9606': 'HLA',
  '9615': 'DLA',
  '9796': 'ELA',
  '9823': 'SLA',
  '9913': 'BoLA',
  '9940': 'Ovar',
  '10090': 'H2',
  '10116': 'rat'
}

def clean_code(name):
  return name.replace('BF-','') \
             .replace('BoLA-','') \
             .replace('Caja-','') \
             .replace('DLA-','') \
             .replace('ELA-','') \
             .replace('Gogo-','') \
             .replace('H2-','') \
             .replace('HLA-','') \
             .replace('Mamu-','') \
             .replace('Papa-','') \
             .replace('Patr-','') \
             .replace('Saoe-','') \
             .replace('SLA-','') \
             .replace('RT1-','')

# Grab the first row and use those headers.

row = csv.reader(args.alleles)
headers = next(row)

# Put the following rows in dicts.

rows = csv.DictReader(args.alleles, fieldnames=headers)

# Update the dicts selectively and print in order.

results = []
for row in rows:
  label = row['displayed_restriction']
  synonyms = row['synonyms'].split(', ')
  includes = list(remove_synonyms(label, synonyms))
  synonyms.sort(key=natural_sort_key)
  includes.sort(key=natural_sort_key)
  row['synonyms'] = '|'.join(synonyms)
  row['includes'] = '|'.join(includes)
  row['organism'] = organisms[row['organism']]
  row['class'] = row['class'].replace('MHC class ', '').replace(' MHC', '').replace('non-','non ')
  row['haplotype'] = row['haplotype'].replace(' haplotype', '')
  row['serotype'] = clean_code(row['serotype'].replace(' serotype', ''))
  row['chain_i_name'] = row['chain_i_name'].replace(' chain', '')
  row['chain_ii_name'] = row['chain_ii_name'].replace(' chain', '')
  row['chain_i_locus'] = clean_code(row['chain_i_locus'].replace(' locus', ''))
  row['chain_ii_locus'] = clean_code(row['chain_ii_locus'].replace(' locus', ''))

  if row['restriction_level'] == 'class':
    row['displayed_restriction'] = codes[row['organism_ncbi_tax_id']] \
        + ' class ' + row['class']

  values = []
  for header in headers:
    values.append(row[header] or '')
  results.append(values)

results.sort(key=lambda x: int(x[0]))

print('\t'.join(headers))
for result in results:
  print('\t'.join(result))

