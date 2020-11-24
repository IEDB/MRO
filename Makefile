# MRO Makefile
# Becky Jackson <rbca.jackson@gmail.com>
#
# WARN: This file contains significant whitespace, i.e. tabs!
# Ensure that your text editor shows you those characters.

### Workflow
#
# Tasks to edit and release MRO.
#
# #### Edit
#
# 1. [Edit Google Sheet](./src/scripts/cogs.sh)
# 2. [Validate Tables](validate_tables) (IDs not required for all terms)
# 3. [Assign New IDs](assign_ids)
# 4. [Validate Tables](validate_tables_strict) (all terms must have IDs)
# 5. [Prepare Products](prepare)
# 6. [Run Tests](test)
# 7. [View Term Table](build/mro.html) or [Browse Tree](./src/scripts/tree.sh)
#
# #### Commit Changes
#
# 1. Run `Status` to see changes
# 2. Run `Commit` and enter message
# 3. Run `Push` and create a new Pull Request
#
# #### Before you go...
#
# [Clean Build Directory](clean)
# [Destroy Google Sheet](destroy)


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
TODAY := $(shell date +%Y-%m-%d)

tables = external core genetic-locus haplotype serotype chain molecule haplotype-molecule serotype-molecule mutant-molecule evidence chain-sequence
source_files = $(foreach o,$(tables),ontology/$(o).tsv)
build_files = $(foreach o,$(tables),build/$(o).tsv)
templates = $(foreach i,$(build_files),--template $(i))


### Set Up

build build/validate:
	mkdir -p $@

# We use the official development version of ROBOT for most things.

build/robot.jar: | build
	curl -L -o $@ https://build.obolibrary.io/job/ontodev/job/robot/job/master/lastSuccessfulBuild/artifact/bin/robot.jar

# Download rdftab based on operating system

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	RDFTAB_URL := https://github.com/ontodev/rdftab.rs/releases/download/v0.1.1/rdftab-x86_64-apple-darwin
else
	RDFTAB_URL := https://github.com/ontodev/rdftab.rs/releases/download/v0.1.1/rdftab-x86_64-unknown-linux-musl
endif

build/rdftab: | build
	curl -L -o $@ $(RDFTAB_URL)
	chmod +x $@


### COGS Tasks

ALL_SHEETS := index.tsv iedb/iedb.tsv $(source_files)
COGS_SHEETS := $(foreach S,$(ALL_SHEETS),.cogs/$(notdir $(S)))

.PHONY: load
load: $(COGS_SHEETS)

.cogs/index.tsv: index.tsv | .cogs
	cogs add $< -r 2

.cogs/iedb.tsv: iedb/iedb.tsv | .cogs
	cogs add $< -r 2

.cogs/%.tsv: ontology/%.tsv | .cogs
	cogs add $< -r 2

.PHONY: push
push: | .cogs
	cogs push

.PHONY: destroy 
destroy: | .cogs
	cogs delete -f


### Table Validation

# Validate the contents of the templates
.PRECIOUS: build/validation_errors.tsv
build/validation_errors.tsv: src/scripts/validate_templates.py index.tsv iedb/iedb.tsv $(build_files)
	python3 $< index.tsv iedb/iedb.tsv build $@ -a

.PRECIOUS: build/validation_errors_strict.tsv
build/validation_errors_strict.tsv: src/scripts/validate_templates.py index.tsv iedb/iedb.tsv $(build_files)
	python3 $< index.tsv iedb/iedb.tsv build $@

apply_%: build/validation_%.tsv | .cogs
	cogs clear all
	cogs apply $<

.PHONY: validate_tables
validate_tables:
	cogs fetch && cogs pull
	make apply_errors
	cogs push

.PHONY: validate_tables_strict
validate_tables_strict:
	cogs fetch && cogs pull
	make apply_errors_strict
	cogs push

### Processing

.PHONY: assign_ids
assign_ids: src/scripts/assign-ids.py index.tsv iedb/iedb.tsv $(source_files)
	# cogs fetch && cogs pull
	python3 $< index.tsv iedb/iedb.tsv ontology
	cogs push

### Review

build/mro.db: src/prefixes.sql mro.owl | build/rdftab
	rm -rf $@
	sqlite3 $@ < $<
	./build/rdftab $@ < $(word 2,$^)

build/mro.html: mro.owl | build/robot.jar
	$(ROBOT) export --input $< --header "ID|LABEL|SubClass Of|definition" --format html --export $@


### Ontology Source Tables

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
	python3 src/sort.py index.tsv $(source_files)

# Check for whitespace during update-tsv step
.PRECIOUS: build/whitespace.tsv
build/whitespace.tsv: src/detect_whitespace.py index.tsv iedb/iedb.tsv iedb/iedb-manual.tsv $(source_files)
	python3 $^ $@


### Sequences

# Moved to GitHub
# OLD: https://github.com/ANHIG/IMGTHLA/raw/Latest/hla_prot.fasta
build/hla.fasta: | build
	curl -o $@ -L https://github.com/ANHIG/IMGTHLA/raw/Latest/hla_prot.fasta

build/mhc.fasta: | build
	wget -O $@ ftp://ftp.ebi.ac.uk/pub/databases/ipd/mhc/MHC_prot.fasta

# update-seqs will only write seqs to terms without seqs
.PHONY: update-seqs
update-seqs: src/update_seqs.py ontology/chain-sequence.tsv build/hla.fasta build/mhc.fasta
	python3 $^

# refresh-seqs will overwrite existing seqs with new seqs
.PHONY: refresh-seqs
refresh-seqs: src/update_seqs.py ontology/chain-sequence.tsv build/hla.fasta build/mhc.fasta
	python3 $^ -o

build/hla_prot.fasta: | build
	curl -o $@ -L https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_prot.fasta

build/AlleleList.txt: | build
	curl -o $@ -L https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt

.PHONY: update-alleles
update-alleles: src/update_human_alleles.py ontology/chain-sequence.tsv ontology/chain.tsv ontology/molecule.tsv ontology/genetic-locus.tsv index.tsv build/hla_prot.fasta build/AlleleList.txt
	python3 $^

.PHONY: update-cow-alleles
update-cow-alleles: src/update_cow_alleles.py ontology/chain-sequence.tsv ontology/chain.tsv ontology/molecule.tsv ontology/genetic-locus.tsv index.tsv build/mhc.fasta iedb/iedb.tsv
	python3 $^

### OWL Files

mro.owl: build/mro-import.owl index.tsv $(build_files) ontology/metadata.ttl | build/robot.jar
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
	--version-iri "$(OBO)/mro/$(TODAY)/mro.owl" \
	--annotation owl:versionInfo "$(TODAY)" \
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
build/mro-iedb.owl: mro.owl iedb/iedb.tsv iedb/iedb-manual.tsv | build/robot.jar iedb
	$(ROBOT) template \
	--prefix "MRO: $(OBO)/MRO_" \
	--input $< \
	--template $(word 2,$^) \
	--template $(word 3,$^) \
	--merge-before \
	--output $@

build/mhc_allele_restriction.csv: build/mro-iedb.owl src/mhc_allele_restriction.rq | build/robot.jar
	$(ROBOT) query --input $(word 1,$^) --select $(word 2,$^) $@

build/mhc_allele_restriction.tsv: src/clean.py build/mhc_allele_restriction.csv ontology/external.tsv | iedb
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

build/mro-base.owl: mro.owl | build/robot.jar
	$(ROBOT) remove --input $< \
	--base-iri $(OBO)/MRO_ \
	--axioms external \
	--output $@

# Run a series of standard OBO checks
.PRECIOUS: build/report.csv
build/report.csv: build/mro-base.owl | build/robot.jar
	$(ROBOT) report --input $< --print 10 --output $@

# Run a series of MRO-specific queries as checks
.PHONY: verify
verify: build/mro-iedb.owl $(VERIFY_QUERIES) | build/robot.jar
	$(ROBOT) verify --input $< \
	--queries $(VERIFY_QUERIES) \
	--output-dir build

# Validate the MHC-allele restriction table
.PRECIOUS: build/mhc_allele_restriction_errors.tsv
build/mhc_allele_restriction_errors.tsv: src/validate_mhc_allele_restriction.py build/mhc_allele_restriction.tsv | build
	python3 $^ $@

.PHONY: test
test: build/validation_errors_strict.tsv
test: build/report.csv
test: verify
test: build/mhc_allele_restriction_errors.tsv

# Python testing
.PHONY: pytest
pytest:
	py.test src/tree.py
	py.test src/synonyms.py


### General

# Prepare products for testing & review
.PHONY: prepare
prepare: build/mro.db
prepare: update-seqs
prepare: update-iedb
prepare:
	pip install -r requirements.txt

.PHONY: clean
clean:
	rm -rf build
	rm -f iedb.zip
	rm -f mro.owl.gz

.PHONY: all
all: clean test


### Release

# IEDB products
iedb.zip: $(IEDB_TARGETS)
	zip -rj $@ $^

# Provide all commits since last tag (excluding merges)
.PHONY: build/release-notes.txt
build/release-notes.txt: | build
	rm -f $@
	echo "New in this release:" >> $@
	git log $$(git describe --tags --abbrev=0)..HEAD --no-merges --oneline \
	| sed "s/^/* /" >> $@

# Release using GitHub CLI
# GITHUB_TOKEN env variable must be set to a PAT with "repo" permissions
# https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/creating-a-personal-access-token
.PHONY: release
release: mro.owl iedb.zip build/release-notes.txt
	gh release create v$(TODAY) mro.owl iedb.zip \
	-t "$(TODAY) Release" \
	-F build/release-notes.txt
