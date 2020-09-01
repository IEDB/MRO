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

build build/validate:
	mkdir -p $@


### ROBOT
#
# We use the official development version of ROBOT for most things.

build/robot.jar: | build
	curl -L -o $@ https://build.obolibrary.io/job/ontodev/job/robot/job/add_validate_operation/lastSuccessfulBuild/artifact/bin/robot.jar


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

update-tsv: update-tsv-files build/whitespace.tsv

# Update TSV files from Excel
.PHONY: update-tsv-files
update-tsv-files:
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

build/mhc.fasta: | build
	wget -O $@ ftp://ftp.ebi.ac.uk/pub/databases/ipd/mhc/MHC_prot.fasta

.PHONY: update-seqs
update-seqs: src/update_seqs.py ontology/chain-sequence.tsv build/hla.fasta build/mhc.fasta
	python3 $^

build/hla_prot.fasta: | build
	curl -o $@ -L https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_prot.fasta

build/AlleleList.txt: | build
	curl -o $@ -L https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt

.PHONY: update-alleles
update-alleles: src/check_missing_alleles.py ontology/chain-sequence.tsv ontology/chain.tsv ontology/molecule.tsv ontology/genetic-locus.tsv index.tsv build/hla_prot.fasta build/AlleleList.txt
	python3 $^

### OWL Files

build/mro.owl: build/mro-import.owl index.tsv $(build_files) ontology/metadata.ttl | build/robot.jar
	$(ROBOT) merge \
	--input $< \
	template \
	--prefix "MRO: $(OBO)/MRO_" \
	--prefix "REO: $(OBO)/REO_" \
	--template index.tsv \
	$(templates) \
	--merge-before \
	reason \
	--reasoner ELK \
	--remove-redundant-subclass-axioms false \
	annotate \
	--ontology-iri "$(OBO)/mro.owl" \
	--version-iri "$(OBO)/mro/$(shell date +%Y-%m-%d)/mro.owl" \
	--annotation owl:versionInfo "$(shell date +%Y-%m-%d)" \
	--annotation-file ontology/metadata.ttl \
	--output $@

build/mro-import.owl: ontology/import.txt $(LIB)/ro.owl $(LIB)/obi.owl $(LIB)/eco.owl | build/robot.jar
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

IEDB_TARGETS := build/mro-iedb.owl \
                build/mhc_allele_restriction.tsv \
                build/ALLELE_FINDER_NAMES.csv \
                build/ALLELE_FINDER_SEARCH.csv \
                build/ALLELE_FINDER_TREE.csv

# extended version for IEDB use
build/mro-iedb.owl: build/mro.owl iedb/iedb.tsv iedb/iedb-manual.tsv | build/robot.jar iedb
	$(ROBOT) template \
	--prefix "MRO: $(OBO)/MRO_" \
	--input $< \
	--template $(word 2,$^) \
	--template $(word 3,$^) \
	--merge-before \
	--output $@

build/mhc_allele_restriction.csv: build/mro-iedb.owl src/mhc_allele_restriction.rq | build/robot.jar
	$(ROBOT) query --input $(word 1,$^) --select $(word 2,$^) $@

build/mhc_allele_restriction.tsv: src/clean.py build/mhc_allele_restriction.csv | iedb
	python3 $^ > $@

build/ALLELE_FINDER_NAMES.csv: build/mro-iedb.owl src/names.rq | build/robot.jar iedb
	$(ROBOT) query --input $(word 1,$^) --select $(word 2,$^) $@.tmp --format csv
	tail -n+2 $@.tmp | dos2unix > $@
	rm $@.tmp

build/ALLELE_FINDER_SEARCH.csv: build/mro-iedb.owl src/search.rq | build/robot.jar iedb
	$(ROBOT) query --input $(word 1,$^) --select $(word 2,$^) $@.tmp --format csv
	tail -n+2 $@.tmp | dos2unix > $@
	rm $@.tmp

build/parents.csv: build/mro-iedb.owl src/parents.rq | build/robot.jar
	$(ROBOT) query --input $(word 1,$^) --select $(word 2,$^) $@

build/ALLELE_FINDER_TREE.csv: src/tree.py build/parents.csv | iedb
	python3 $^ --mode CSV > $@

build/tree.json: src/tree.py build/parents.csv | build
	python3 $^ --mode JSON > $@

build/full_tree.json: src/tree.py build/full_tree.csv | build
	python3 $^ --mode JSON > $@

.PHONY: update-iedb
update-iedb: $(IEDB_TARGETS)


### Testing & verification

VERIFY_QUERIES = $(wildcard src/verify/*.rq)

build/mro-base.owl: build/mro.owl | build/robot.jar
	$(ROBOT) remove --input $< \
	--base-iri $(OBO)/MRO_ \
	--axioms external \
	--output $@

.PRECIOUS: build/report.csv
build/report.csv: build/mro-base.owl | build/robot.jar
	$(ROBOT) report --input $< --print 10 --output $@

.PHONY: verify
verify: build/mro-iedb.owl $(VERIFY_QUERIES) | build/robot.jar
	$(ROBOT) verify --input $< \
	--queries $(VERIFY_QUERIES) \
	--output-dir build

# Validate a relaxed/reduced version of MRO
.PHONY: validate
validate: build/mro.owl $(source_files) | build/robot.jar build/validate
	$(ROBOT) relax --input $< \
	reduce remove --axioms equivalent \
	validate $(foreach i,$(source_files),--table $(i)) \
	--skip-row 2 \
	--format html \
	--output-dir build/validate

.PRECIOUS: build/mhc_allele_restriction_errors.tsv
build/mhc_allele_restriction_errors.tsv: src/validate_mhc_allele_restriction.py build/mhc_allele_restriction.tsv | build
	python3 $^ $@

.PRECIOUS: build/whitespace.tsv
build/whitespace.tsv: src/detect_whitespace.py index.tsv iedb/iedb.tsv iedb/iedb-manual.tsv $(source_files)
	python3 $^ $@


### Release files

iedb.zip: $(IEDB_TARGETS)
	zip $@ $^

mro.owl.gz: build/mro.owl
	gzip -f $<

release: iedb.zip mro.owl.gz


### General

.PHONY: test
test: build/report.csv verify build/mhc_allele_restriction_errors.tsv

.PHONY: pytest
pytest:
	py.test src/tree.py
	py.test src/synonyms.py

.PHONY: clean
clean:
	rm -rf build
	rm -f iedb.zip
	rm -f mro.owl.gz

.PHONY: all
all: clean release
