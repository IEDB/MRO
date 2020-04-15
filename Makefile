### MHC Restriction Ontology Makefile
#
# James A. Overton <james@overton.ca>
#
# This file is used to build MRO from source.
# Usually you want to run:
#
#     make clean all
#
# Requirements:
#
# - GNU Make
# - ROBOT <http://github.com/ontodev/robot>
# - Python 3
#   - openpyxl module <https://openpyxl.readthedocs.io>


### Configuration
#
# These are standard options to make Make sane:
# <http://clarkgrubb.com/makefile-style-guide#toc2>

MAKEFLAGS += --warn-undefined-variables
SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
.SUFFIXES:
.SECONDARY:

OBO = http://purl.obolibrary.org/obo
LIB = lib
ROBOT := java -jar build/robot.jar


### Set Up

build:
	mkdir -p $@


### ROBOT
#
# We use the official development version of ROBOT for most things.

build/robot.jar: | build
	curl -L -o $@ https://github.com/ontodev/robot/releases/download/v1.6.0/robot.jar


### Ontology Source Tables

tables = external core genetic-locus haplotype serotype chain molecule haplotype-molecule serotype-molecule mutant-molecule evidence chain-sequence
source_files = $(foreach o,$(tables),ontology/$(o).tsv)
build_files = $(foreach o,$(tables),build/$(o).tsv)
templates = $(foreach i,$(build_files),--template $(i))

build/%.tsv: ontology/%.tsv | build
	cp $< $@

# Generate automatic synonyms
build/molecule.tsv: src/synonyms.py ontology/molecule.tsv | build
	python3 $^ > $@
build/haplotype-molecule.tsv: src/synonyms.py ontology/haplotype-molecule.tsv | build
	python3 $^ > $@
build/serotype-molecule.tsv: src/synonyms.py ontology/serotype-molecule.tsv | build
	python3 $^ > $@
build/mutant-molecule.tsv: src/synonyms.py ontology/mutant-molecule.tsv | build
	python3 $^ > $@

# Represent tables in Excel
mro.xlsx: src/tsv2xlsx.py index.tsv iedb/iedb.tsv ontology/genetic-locus.tsv ontology/haplotype.tsv ontology/serotype.tsv ontology/chain.tsv ontology/chain-sequence.tsv ontology/molecule.tsv ontology/haplotype-molecule.tsv ontology/serotype-molecule.tsv ontology/mutant-molecule.tsv ontology/core.tsv ontology/external.tsv iedb/iedb-manual.tsv ontology/evidence.tsv
	python3 $< $@ $(wordlist 2,100,$^)

# Update TSV files from Excel
.PHONY: update-tsv
update-tsv:
	python3 src/xlsx2tsv.py mro.xlsx index > index.tsv
	python3 src/xlsx2tsv.py mro.xlsx iedb > iedb/iedb.tsv
	python3 src/xlsx2tsv.py mro.xlsx iedb-manual > iedb/iedb-manual.tsv
	$(foreach t,$(tables),python3 src/xlsx2tsv.py mro.xlsx $(t) > ontology/$(t).tsv;)
	python3 src/sort.py $(source_files)

# Sort TSV files by first column
.PHONY: sort
sort:
	python3 src/sort.py $(source_files)


### Sequences

# Moved to GitHub
# OLD: https://github.com/ANHIG/IMGTHLA/raw/Latest/hla_prot.fasta
build/hla.fasta: | build
	curl -o $@ -L https://github.com/ANHIG/IMGTHLA/raw/Latest/hla_prot.fasta

build/sla.fasta: | build
	curl -o $@ -L https://www.ebi.ac.uk/ipd/mhc/group/SLA/download?type=protein

build/dla.fasta: | build
	curl -o $@ -L https://www.ebi.ac.uk/ipd/mhc/group/DLA/download?type=protein

.PHONY: update-seqs
update-seqs: src/update_seqs.py ontology/chain-sequence.tsv build/hla.fasta build/sla.fasta build/dla.fasta
	python3 $^


### OWL Files

mro.owl: mro-import.owl index.tsv $(build_files) ontology/metadata.ttl | build/robot.jar
	$(ROBOT) merge \
	--input mro-import.owl \
	template \
	--prefix "MRO: $(OBO)/MRO_" \
	--prefix "REO: $(OBO)/REO_" \
	--template index.tsv \
	$(templates) \
	--merge-before \
	reason --reasoner HermiT \
	--remove-redundant-subclass-axioms false \
	annotate \
	--ontology-iri "$(OBO)/mro.owl" \
	--version-iri "$(OBO)/mro/$(shell date +%Y-%m-%d)/mro.owl" \
	--annotation owl:versionInfo "$(shell date +%Y-%m-%d)" \
	--annotation-file ontology/metadata.ttl \
	--output $@

mro-import.owl: ontology/import.txt $(LIB)/ro.owl $(LIB)/obi.owl $(LIB)/eco.owl | build/robot.jar
	$(ROBOT) merge \
	--input $(LIB)/eco.owl \
	--input $(LIB)/obi.owl \
	--input $(LIB)/ro.owl \
	extract \
	--method MIREOT \
	--prefix "REO: $(OBO)/REO_" \
	--upper-term "GO:0008150" \
	--upper-term "IAO:0000030" \
	--upper-term "OBI:1110128" \
	--upper-term "ECO:0000000" \
	--upper-term "BFO:0000040" \
	--upper-term "PR:000000001" \
	--lower-terms $< \
	--output $@

# fetch ontology dependencies
$(LIB)/%:
	mkdir -p $(LIB)
	cd $(LIB) && curl -LO "$(OBO)/$*"


### Generate files for IEDB
#
# This includes an extended OWL file
# and tables for the IEDB database and Finders.

iedb:
	mkdir -p $@

# extended version for IEDB use
iedb/mro-iedb.owl: mro.owl iedb/iedb.tsv iedb/iedb-manual.tsv | build/robot.jar iedb
	$(ROBOT) template \
	--prefix "MRO: $(OBO)/MRO_" \
	--input $< \
	--template $(word 2,$^) \
	--template $(word 3,$^) \
	--merge-before \
	--output $@

build/mhc_allele_restriction.csv: iedb/mro-iedb.owl src/mhc_allele_restriction.rq | build/robot.jar
	$(ROBOT) query --input $(word 1,$^) --select $(word 2,$^) $@

iedb/mhc_allele_restriction.tsv: src/clean.py build/mhc_allele_restriction.csv | iedb
	python3 $^ > $@

iedb/ALLELE_FINDER_NAMES.csv: iedb/mro-iedb.owl src/names.rq | build/robot.jar iedb
	$(ROBOT) query --input $(word 1,$^) --select $(word 2,$^) $@.tmp --format csv
	tail -n+2 $@.tmp | dos2unix > $@
	rm $@.tmp

iedb/ALLELE_FINDER_SEARCH.csv: iedb/mro-iedb.owl src/search.rq | build/robot.jar iedb
	$(ROBOT) query --input $(word 1,$^) --select $(word 2,$^) $@.tmp --format csv
	tail -n+2 $@.tmp | dos2unix > $@
	rm $@.tmp

build/parents.csv: iedb/mro-iedb.owl src/parents.rq | build/robot.jar
	$(ROBOT) query --input $(word 1,$^) --select $(word 2,$^) $@

iedb/ALLELE_FINDER_TREE.csv: src/tree.py build/parents.csv | iedb
	python3 $^ --mode CSV > $@

build/tree.json: src/tree.py build/parents.csv | build
	python3 $^ --mode JSON > $@

build/full_tree.json: src/tree.py build/full_tree.csv | build
	python3 $^ --mode JSON > $@

IEDB_TARGETS := iedb/mro-iedb.owl iedb/mhc_allele_restriction.tsv iedb/ALLELE_FINDER_NAMES.csv iedb/ALLELE_FINDER_SEARCH.csv iedb/ALLELE_FINDER_TREE.csv
.PHONY: update-iedb
update-iedb: $(IEDB_TARGETS)


### General

VERIFY_QUERIES = $(wildcard src/verify/*.rq)

build/mro-base.owl: mro.owl | build/robot.jar
	$(ROBOT) remove --input $< \
	--base-iri $(OBO)/MRO_ \
	--axioms external \
	--output $@

.PRECIOUS: build/report.csv
build/report.csv: build/mro-base.owl | build/robot.jar
	$(ROBOT) report --input $< \
	--output $@

.PHONY: verify
verify: iedb/mro-iedb.owl $(VERIFY_QUERIES) | build/robot.jar
	$(ROBOT) verify --input $< \
	--queries $(VERIFY_QUERIES) \
	--output-dir build

.PHONY: test
test: build/report.csv verify

.PHONY: pytest
pytest:
	py.test src/tree.py
	py.test src/synonyms.py

.PHONY: clean
clean:
	rm -rf build
	rm -f mro.owl mro-import.owl
	rm -f $(IEDB_TARGETS)

.PHONY: all
all: update-iedb
