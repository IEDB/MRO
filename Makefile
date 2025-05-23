# MRO Makefile
# Becky Jackson <rbca.jackson@gmail.com>
#
# WARN: This file contains significant whitespace, i.e. tabs!
# Ensure that your text editor shows you those characters.

### Workflow
#
# Tasks to edit and release MRO.
#
# * [Upload mro.xlsx](./src/scripts/sheet.py?action=create)
# * [Download mro.xlsx](mro.xlsx)
#
# #### Build Products
#
# 1. [Validate Tables](validate_tables) (IDs not required for all terms)
# 2. [Assign New IDs](assign_ids)
# 3. [Validate Tables](validate_tables_strict) (all terms must have IDs)
# 4. [Prepare Products](prepare)
# 5. View the results:
#     * [Table Diffs](build/diff.html)
#     * [Term Table](build/mro.html)
#     * [Browse Tree](./src/scripts/tree.sh)
# 6. [Run Tests](test)
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

build build/validate build/master build/diff:
	mkdir -p $@

# We use the official development version of ROBOT for most things.

build/robot.jar: | build
	curl -L -o $@ https://github.com/ontodev/robot/releases/download/v1.8.1/robot.jar

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


### Table Validation

# Validate the contents of the templates
.PRECIOUS: build/validation_errors.tsv
build/validation_errors.tsv: src/scripts/validation/validate_templates.py index.tsv iedb/iedb.tsv $(source_files) | build
	python3 $< index.tsv iedb/iedb.tsv ontology $@ -a

.PRECIOUS: build/validation_errors_strict.tsv
build/validation_errors_strict.tsv: src/scripts/validation/validate_templates.py index.tsv iedb/iedb.tsv $(source_files) | build
	python3 $< index.tsv iedb/iedb.tsv ontology $@ -f

apply_%: build/validation_%.tsv
	axle clear all
	axle apply $<

.PHONY: validate_tables
validate_tables:
	make apply_errors
	axle push

.PHONY: validate_tables_strict
validate_tables_strict:
	make apply_errors_strict
	axle push

### Processing

.PHONY: assign_ids
assign_ids: src/scripts/assign-ids.py index.tsv iedb/iedb.tsv $(source_files)
	python3 $< index.tsv iedb/iedb.tsv ontology
	axle push

### Review

build/mro.db: src/queries/prefixes.sql mro.owl | build/rdftab
	rm -rf $@
	sqlite3 $@ < $<
	./build/rdftab $@ < $(word 2,$^)
	sqlite3 $@ "CREATE INDEX idx_stanza ON statements (stanza);"
	sqlite3 $@ "CREATE INDEX idx_subject ON statements (subject);"
	sqlite3 $@ "CREATE INDEX idx_predicate ON statements (predicate);"
	sqlite3 $@ "CREATE INDEX idx_object ON statements (object);"
	sqlite3 $@ "CREATE INDEX idx_value ON statements (value);"
	sqlite3 $@ "ANALYZE;"

build/mro.html: mro.owl | build/robot.jar
	$(ROBOT) export --input $< --header "ID|LABEL|SubClass Of|definition" --format html --export $@


### Ontology Source Tables

# Replace labels for any term with single quotes used in a C ROBOT template string
build/%-fixed.tsv: src/scripts/replace_labels.py ontology/%.tsv | build
	python3 $^ $@

# Some templates don't need automatic synonyms, just copy these directly
build/%.tsv: build/%-fixed.tsv | build
	cp $< $@

# Generate automatic synonyms for these tables
build/molecule.tsv: src/scripts/synonyms.py build/molecule-fixed.tsv | build
	python3 $^ > $@
build/haplotype-molecule.tsv: src/scripts/synonyms.py build/haplotype-molecule-fixed.tsv | build
	python3 $^ > $@
build/serotype-molecule.tsv: src/scripts/synonyms.py build/serotype-molecule-fixed.tsv | build
	python3 $^ > $@
build/mutant-molecule.tsv: src/scripts/synonyms.py build/mutant-molecule-fixed.tsv | build
	python3 $^ > $@

# Represent tables in Excel
mro.xlsx: src/scripts/tsv2xlsx.py index.tsv iedb/iedb.tsv ontology/genetic-locus.tsv ontology/haplotype.tsv ontology/serotype.tsv ontology/chain.tsv ontology/chain-sequence.tsv ontology/molecule.tsv ontology/haplotype-molecule.tsv ontology/serotype-molecule.tsv ontology/mutant-molecule.tsv ontology/core.tsv ontology/external.tsv iedb/iedb-manual.tsv ontology/evidence.tsv ontology/rejected.tsv
	python3 $< $@ $(wordlist 2,100,$^)

update-tsv: update-tsv-files sort build/whitespace.tsv

# Update TSV files from Excel
.PHONY: update-tsv-files
update-tsv-files:
	python3 src/scripts/xlsx2tsv.py mro.xlsx index > index.tsv
	python3 src/scripts/xlsx2tsv.py mro.xlsx iedb > iedb/iedb.tsv
	python3 src/scripts/xlsx2tsv.py mro.xlsx iedb-manual > iedb/iedb-manual.tsv
	$(foreach t,$(tables) rejected,python3 src/scripts/xlsx2tsv.py mro.xlsx $(t) > ontology/$(t).tsv;)
	python3 src/scripts/sort.py $(source_files)

# Sort TSVs
.PHONY: sort
sort: sort-templates sort-iedb

# Sort template files by first column
.PHONY: sort-templates
sort-templates:
	python3 src/scripts/sort.py index.tsv $(source_files)

# Sort IEDB table by IEDB ID
.PHONY: sort-iedb
sort-iedb: src/scripts/sort_iedb.py iedb/iedb.tsv
	python3 $^

# Check for whitespace during update-tsv step
.PRECIOUS: build/whitespace.tsv
build/whitespace.tsv: src/scripts/validation/detect_whitespace.py index.tsv iedb/iedb.tsv iedb/iedb-manual.tsv $(source_files)
	python3 $^ $@


### Sequences

# Moved to GitHub
# OLD: https://github.com/ANHIG/IMGTHLA/raw/Latest/hla_prot.fasta
build/hla.fasta: | build
	curl -L -o $@ https://github.com/ANHIG/IMGTHLA/raw/Latest/hla_prot.fasta

build/mhc.fasta: | build
	curl -L -o $@ ftp://ftp.ebi.ac.uk/pub/databases/ipd/mhc/MHC_prot.fasta

# update-seqs will only write seqs to terms without seqs
.PHONY: update-seqs
update-seqs: src/scripts/update_seqs.py ontology/chain-sequence.tsv build/hla.fasta build/mhc.fasta
	python3 $^

# refresh-seqs will overwrite existing seqs with new seqs
.PHONY: refresh-seqs
refresh-seqs: src/scripts/update_seqs.py ontology/chain-sequence.tsv build/hla.fasta build/mhc.fasta
	python3 $^ -o

# refresh-hla-seqs will overwrite existing HLA seqs with new seqs
.PHONY: refresh-hla-seqs
refresh-hla-seqs: src/scripts/update_seqs.py ontology/chain-sequence.tsv build/hla.fasta
	python3 $^ -o -H

# Diff chain sequences against `master` as a table.
build/chain-sequence-diff.tsv: ontology/chain-sequence.tsv
	head -n1 $< > $@
	git diff --unified=0 master | grep "^[+-]HLA-" >> $@

build/hla_prot.fasta: | build
	curl -o $@ -L https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_prot.fasta

build/AlleleList.txt: | build
	curl -o $@ -L https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt

.PHONY: update-mhcflurry-alleles
update-mhcflurry-alleles: src/scripts/alleles/update_alleles_mhcflurry.py ontology/chain-sequence.tsv ontology/chain.tsv ontology/molecule.tsv ontology/genetic-locus.tsv index.tsv build/mhc.fasta iedb/iedb.tsv
	python3 $^

### OWL Files

mro.owl: build/mro-import.owl index.tsv $(build_files) ontology/metadata.ttl | build/robot.jar
	$(ROBOT) template \
	--input $< \
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

build/mro-import.owl: build/eco-import.ttl build/iao-import.ttl build/obi-import.ttl build/ro-import.ttl ontology/import.txt | build/robot.jar
	$(ROBOT) merge \
	--input build/eco-import.ttl \
	--input build/obi-import.ttl \
	--input build/ro-import.ttl \
	--input build/iao-import.ttl \
	extract \
	--method MIREOT \
	--upper-term "GO:0008150" \
	--upper-term "IAO:0000030" \
	--upper-term "OBI:1110128" \
	--upper-term "ECO:0000000" \
	--upper-term "BFO:0000040" \
	--upper-term "PR:000000001" \
	--lower-terms $(word 5,$^) \
	--output $@

# fetch ontology dependencies
$(LIB)/%:
	mkdir -p $(LIB)
	cd $(LIB) && curl -LO "$(OBO)/$*"

UC = $(shell echo '$1' | tr '[:lower:]' '[:upper:]')

# OBI IAO:0000115 has mulitples so get the definiton from here
# we could also just add this to index.tsv
build/%.txt: ontology/import.txt | build
	sed -n '/$(call UC,$(notdir $(basename $@)))/p' $< > $@

# RO:0000056 isn't in RO?
# we could also just add this to index.tsv
build/obi.txt: ontology/import.txt | build
	sed '/^ECO/d' $< | sed '/^RO/d' | sed '/^IAO/d' > $@
	echo "RO:0000056" >> $@

build/%.db: src/queries/prefixes.sql $(LIB)/%.owl | build/rdftab
	rm -rf $@
	sqlite3 $@ < $<
	./build/rdftab $@ < $(word 2,$^)

PREDICATES := IAO:0000111 IAO:0000112 IAO:0000115 IAO:0000119 IAO:0000412 OBI:9991118 rdfs:label

build/%-import.ttl: build/%.db build/%.txt
	python3 -m gizmos.extract -d $< -T $(word 2,$^) $(foreach P,$(PREDICATES), -p $(P)) > $@

### Generate files for IEDB
#
# This includes an extended OWL file
# and tables for the IEDB database and Finders.

IEDB_TARGETS := build/mro-iedb.owl \
                build/mhc_allele_restriction.tsv \
                build/molecule_export.tsv \
                build/ALLELE_FINDER_NAMES.csv \
                build/ALLELE_FINDER_SEARCH.csv \
                build/ALLELE_FINDER_TREE.csv

# add chain_i and chain_ii accessions to IEDB sheet
build/iedb.tsv: src/scripts/add_chain_accessions.py index.tsv iedb/iedb.tsv ontology/molecule.tsv ontology/chain-sequence.tsv
	python3 $^ $@

# extended version for IEDB use
build/mro-iedb.owl: mro.owl build/iedb.tsv iedb/iedb-manual.tsv | build/robot.jar iedb
	$(ROBOT) template \
	--prefix "MRO: $(OBO)/MRO_" \
	--input $< \
	--template $(word 2,$^) \
	--template $(word 3,$^) \
	--merge-before \
	--output $@

build/mro-tdb: build/mro-iedb.owl
	rm -rf $@
	$(ROBOT) query --input $< \
	--create-tdb true \
	--tdb-directory $@

build/mhc_allele_restriction.csv: build/mro-tdb src/queries/mhc_allele_restriction.rq | build/robot.jar
	$(ROBOT) query \
	--tdb-directory $< \
	--keep-tdb-mappings true \
	--select $(word 2,$^) $@

build/mhc_allele_restriction.tsv: src/scripts/clean.py build/mhc_allele_restriction.csv ontology/external.tsv | iedb
	python3 $^ > $@

build/molecule_export.tsv: src/scripts/export_molecule.py index.tsv ontology/external.tsv ontology/molecule.tsv
	python3 $^ $@

build/ALLELE_FINDER_NAMES.csv: build/mro-tdb src/queries/names.rq | build/robot.jar iedb
	$(ROBOT) query \
	--tdb-directory $< \
	--keep-tdb-mappings true \
	--select $(word 2,$^) $@.tmp \
	--format csv
	tail -n+2 $@.tmp | dos2unix > $@
	rm $@.tmp

build/ALLELE_FINDER_SEARCH.csv: build/mro-tdb src/queries/search.rq | build/robot.jar iedb
	$(ROBOT) query \
	--tdb-directory $< \
	--keep-tdb-mappings true \
	--select $(word 2,$^) $@.tmp \
	--format csv
	tail -n+2 $@.tmp | dos2unix > $@
	rm $@.tmp

build/parents.csv: build/mro-tdb src/queries/parents.rq | build/robot.jar
	$(ROBOT) query \
	--tdb-directory $< \
	--keep-tdb-mappings true \
	--select $(word 2,$^) $@

build/ALLELE_FINDER_TREE.csv: src/scripts/tree.py build/parents.csv | iedb
	python3 $^ --mode CSV > $@

build/tree.json: src/scripts/tree.py build/parents.csv | build
	python3 $^ --mode JSON > $@

build/full_tree.json: src/scripts/tree.py build/full_tree.csv | build
	python3 $^ --mode JSON > $@

.PHONY: update-iedb
update-iedb: $(IEDB_TARGETS)


### Testing & verification

VERIFY_QUERIES = $(wildcard src/queries/verify/*.rq)

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
build/mhc_allele_restriction_errors.tsv: src/scripts/validation/validate_mhc_allele_restriction.py build/mhc_allele_restriction.tsv | build
	python3 $^ $@

.PHONY: test
test: build/validation_errors_strict.tsv
test: build/report.csv
test: verify
test: build/mhc_allele_restriction_errors.tsv

# Python testing
.PHONY: pytest
pytest:
	py.test src/scripts/tree.py
	py.test src/scripts/synonyms.py

# Table diffs

DIFF_TABLES := build/diff/index.html build/diff/iedb.html $(foreach S,$(tables),build/diff/$(S).html) build/diff/iedb-manual.html

# Workaround to make sure master branch exists
build/fetched.txt:
	git fetch origin master:master && date > $@

build/diff/%.html: ontology/%.tsv build/fetched.txt | build/master build/diff
	git show master:$< > build/master/$(notdir $<)
	daff build/master/$(notdir $<) $< --output $@ --fragment

build/diff/iedb.html: iedb/iedb.tsv build/fetched.txt | build/master build/diff
	git show master:$< > build/master/$(notdir $<)
	daff build/master/$(notdir $<) $< --output $@ --fragment

build/diff/iedb-manual.html: iedb/iedb-manual.tsv build/fetched.txt | build/master build/diff
	git show master:$< > build/master/$(notdir $<)
	daff build/master/$(notdir $<) $< --output $@ --fragment

build/diff/index.html: index.tsv build/fetched.txt | build/master build/diff
	git show master:$< > build/master/$(notdir $<)
	daff build/master/$(notdir $<) $< --output $@ --fragment

build/diff.html: src/scripts/diff.py src/scripts/diff.html $(DIFF_TABLES)
	python3 $< src/scripts/diff.html $(tables) > $@


### General

# Prepare products for testing & review
.PHONY: prepare
prepare: sort
prepare: clean
# prepare: build/mro.db
prepare: update-seqs
prepare: test
prepare: update-iedb
# prepare: build/diff.html

.PHONY: clean
clean:
	rm -rf mro.owl
	rm -rf build
	rm -rf lib
	rm -f iedb.zip
	rm -f mro.owl.gz

.PHONY: all
all: clean prepare


### Release

# IEDB products
iedb.zip: $(IEDB_TARGETS)
	zip -rj $@ $^

# Release using GitHub CLI
# GITHUB_TOKEN env variable must be set to a PAT with "repo" permissions
# https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/creating-a-personal-access-token
.PHONY: release
release: mro.owl mro.xlsx iedb.zip
	# Make sure we're not ahead or behind origin/master
	git branch --show-current | grep master
	git fetch
	! git status -uno | grep -E "(ahead|behind)"
	# Make the release
	gh release create v$(TODAY) \
	--title "$(TODAY) Release" \
	--generate-notes \
	$^
