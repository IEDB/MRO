#!/usr/bin/env python3
#
# Generate a tree representation
# from SPARQL result listing ?subject and ?parent.
# The CSV representation is for IEDB Finder TREE tables.
# The JSON representation works with inspire-tree.

import argparse, csv, re, sys
from collections import defaultdict
from io import StringIO


### Read Rows
#
# All the output methods start with a SPARQL query result table.
# It's usually in CSV format, and we use csv.DictReader.
# No order should be assumed in the SPARQL results!
#
# We read the data and return:
# - data: a dict of dicts
# - roots: a set of IRIs that are used as parents
#   but are not also subjects (children).

example_rows = [
  {'subject': 'd', 'parent': 'c', 'label': 'D'},
  {'subject': 'b', 'parent': 'a', 'label': 'B'},
  {'subject': 'c', 'parent': 'f', 'label': 'C'},
  {'subject': 'f', 'parent': 'a', 'label': 'F'},
  {'subject': 'c', 'parent': 'b', 'label': 'C'},
  {'subject': 'e', 'parent': 'c', 'label': 'E'},
  {'subject': 'a', 'parent': '',  'label': 'A'}
]
example_data = {
  'a': {'label': 'A', 'parents': set(),     'children': {'b', 'f'}},
  'b': {'label': 'B', 'parents': {'a'},     'children': {'c'}},
  'c': {'label': 'C', 'parents': {'b','f'}, 'children': {'d', 'e'}},
  'd': {'label': 'D', 'parents': {'c'},     'children': set()},
  'e': {'label': 'E', 'parents': {'c'},     'children': set()},
  'f': {'label': 'F', 'parents': {'a'},     'children': {'c'}},
}
example_roots = {'a'}

def init_data():
  return {
    'parents': set(),
    'children': set()
  }

def read_rows(rows):
  """Given an iterator of dictionaries,
  return the tuple of the data and roots."""
  data = defaultdict(init_data)
  roots = set()
  for row in rows:
    subject = row['subject']
    parent = row['parent']

    if parent == 'http://www.w3.org/2002/07/owl#Thing':
      parent = ''
    if parent == '':
      roots.add(subject)
    else:
      data[subject]['parents'].add(parent)
      data[parent]['children'].add(subject)

    for key in ['label', 'sort', 'synonyms']:
      if key in row:
        data[subject][key] = row[key]

    if 'synonyms' in data[subject]:
      synonyms = data[subject]['synonyms'].split(', ')
      sorted(synonyms)
      data[subject]['synonyms'] = ', '.join(synonyms)

  return data, roots

def test_read_rows():
  data, roots = read_rows(example_rows)
  assert example_data == data
  assert example_roots == roots

def get_label(data, iri):
  return data[iri]['label']

def get_sort_label(data, iri):
  if 'sort' in data[iri] and data[iri]['sort'] != '':
    return data[iri]['sort']
  else:
    return data[iri]['label']

def get_synonyms(data, iri):
  if 'synonyms' in data[iri]:
    return data[iri]['synonyms']
  else:
    return ''

# http://stackoverflow.com/a/16090640
def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
  return [
    int(text) if text.isdigit() else text.lower()
    for text in re.split(_nsre, s)
  ]

def test_natural_sort():
  natural = ['a1', 'a50', 'a100', 'b200']
  alpha = sorted(natural[:])
  assert natural == sorted(alpha, key=natural_sort_key)

def iris_by_label(data, iris):
  """Given the data dict and a set of IRIs,
  sort the IRIs by their labels using natural numeric order."""
  labels = {}
  for iri in iris:
    labels[get_sort_label(data, iri)] = iri
  result = []
  for label in sorted(labels.keys(), key=natural_sort_key):
    result.append(labels[label])
  return result

def children_by_label(data, parent_iri):
  return iris_by_label(data, data[parent_iri]['children'])


### Write Text
#
# The simplest output is text,
# where each node is indented by a number of dashes.

example_text = """- A
-- B
--- C
---- D
---- E
-- F
--- C
---- D
---- E
"""

def write_text(writer, data, roots):
  for root in iris_by_label(data, roots):
    write_lines(writer, data, root, 1)

def write_lines(writer, data, iri, depth):
  """Recursively write lines starting with dashes
  representing a tree."""
  writer.write('-' * depth + ' ' + get_label(data, iri) + '\n')
  for child in children_by_label(data, iri):
    write_lines(writer, data, child, depth+1)

def test_write_lines():
  output = StringIO()
  write_text(output, example_data, example_roots)
  assert example_text == output.getvalue()


### Write Table
#
# Write a table representation of the tree.
# This is used for IEDB Finders.
# Each row represents a node in the tree,
# and has an index and a parent index.
# If an ontology term has multiple parents,
# it will have multiple rows with distinct indices.

example_csv = """1,a,0,A,,1
2,b,1,B,,1:2
3,c,2,C,,1:2:3
4,d,3,D,,1:2:3:4
5,e,3,E,,1:2:3:5
6,f,1,F,,1:6
7,c,6,C,,1:6:7
8,d,7,D,,1:6:7:8
9,e,7,E,,1:6:7:9
"""

# This global index is used to increment entry numbers.
index = 1

def write_table(writer, data, roots):
  global index
  index = 1
  for root in iris_by_label(data, roots):
    write_rows(writer, data, root, [0])

def write_rows(writer, data, iri, ancestry):
  """Recursively write rows for a 'tree table'."""
  global index
  ancestry = ancestry + [index]
  writer.writerow([
    index,
    iri,
    ancestry[-2],
    get_label(data, iri),
    get_synonyms(data, iri),
    ':'.join(str(x) for x in ancestry[1:])
  ])
  for child in children_by_label(data, iri):
    index += 1
    write_rows(writer, data, child, ancestry)

def test_write_csv():
  output = StringIO()
  writer = csv.writer(output, lineterminator='\n')
  write_table(writer, example_data, example_roots)
  assert example_csv == output.getvalue()


### Write JSON
#
# Output a JSON representation,
# suitable for use with inspire-tree.
# We use a crude string-based method for generating JSON
# to avoid putting it all back in memory.

example_json = """[
 {"text": "A",
  "iri": "a",
  "children": [
  {"text": "B",
   "iri": "b",
   "children": [
   {"text": "C",
    "iri": "c",
    "children": [
    {"text": "D",
     "iri": "d"},
    {"text": "E",
     "iri": "e"}
   ]}
  ]},
  {"text": "F",
   "iri": "f",
   "children": [
   {"text": "C",
    "iri": "c",
    "children": [
    {"text": "D",
     "iri": "d"},
    {"text": "E",
     "iri": "e"}
   ]}
  ]}
 ]}
]
"""

def write_json(writer, data, roots):
  writer.write('[\n')
  sorted_roots = iris_by_label(data, roots)
  for root in sorted_roots:
    write_json_objects(writer, data, root, 1)
    if root != sorted_roots[-1]:
      writer.write(',\n')
  writer.write('\n]\n')

def write_json_objects(writer, data, iri, depth):
  """Recursively write nested JSON objects
  representing a tree."""
  writer.write(' ' * depth + '{"text": "%s",\n' % get_label(data, iri))
  if data[iri]['children']: # node has children
    writer.write(' ' * depth + ' "iri": "%s",\n' % iri)
    writer.write(' ' * depth + ' "children": [\n')
    children = children_by_label(data, iri)
    for child in children:
      write_json_objects(writer, data, child, depth+1)
      if child != children[-1]:
        writer.write(',\n')
      else:
        writer.write('\n')
    writer.write(' ' * depth + ']}')
  else:
    writer.write(' ' * depth + ' "iri": "%s"}' % iri)

def test_write_json():
  output = StringIO()
  write_json(output, example_data, example_roots)
  assert example_json == output.getvalue()


### Main
#
# Parse arguments and run the tool.

def main():
  parser = argparse.ArgumentParser(
      description='Generate Turtle for pruned taxa')
  parser.add_argument('--mode',
      type=str,
      choices=['TEXT','CSV','TSV','JSON'],
      default='TEXT',
      help='output mode')
  parser.add_argument('query',
      type=argparse.FileType('r'),
      help='read query result CSV')
  args = parser.parse_args()

  rows = csv.DictReader(args.query)
  data, roots = read_rows(rows)

  if args.mode == 'TEXT':
    write_text(sys.stdout, data, roots)

  elif args.mode == 'CSV':
    writer = csv.writer(sys.stdout,
        lineterminator='\n')
    write_table(writer, data, roots)

  elif args.mode == 'TSV':
    writer = csv.writer(sys.stdout,
        delimiter='\t',
        lineterminator='\n')
    write_table(writer, data, roots)

  elif args.mode == 'JSON':
    write_json(sys.stdout, data, roots)

if __name__ == "__main__":
  main()


