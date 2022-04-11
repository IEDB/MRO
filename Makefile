# MRO Makefile
# Becky Jackson <rbca.jackson@gmail.com>
#
# WARN: This file contains significant whitespace, i.e. tabs!
# Ensure that your text editor shows you those characters.

### Workflow
#
# 1. [View prototype](./src/scripts/run.py)
# 2. [Update templates](refresh_templates)
# 3. [Rebuild MRO](all)
# 4. [Reload MRO to database](load_mro)
#
# Danger zone:
# * [Refresh database](refresh_db) (WARNING: this takes a long time!)
#

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
ROBOT := java -jar build/robot.jar --prefix "iedb: http://iedb.org/" --prefix "REO: $(OBO)/REO_"
LDTAB := java -jar build/ldtab.jar
TODAY := $(shell date +%Y-%m-%d)

tables = external core genetic-locus haplotype serotype chain molecule haplotype-molecule serotype-molecule mutant-molecule evidence obsolete
source_files = $(foreach o,$(tables),ontology/$(o).tsv)
build_files = $(foreach o,$(tables),build/$(o).tsv)

define \n


endef


### Set Up

build build/validate build/master build/diff build/tables:
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

build/ldtab.jar: | build
	curl -L -o $@ "https://github.com/ontodev/ldtab.clj/releases/download/v2022-03-17/ldtab.jar"


### Ontology Source Tables

# Index containing ID, Label, and Type for all MRO terms & imports
build/index.tsv: $(source_files) | build
	echo -e "ID\tLabel\tType\ttable\nID\tLABEL\tTYPE\t" >> $@
	$(foreach f,$(source_files),tail -n +3 $(f) | cut -f1,2,3 | sed -e 's/$$/\t$(notdir $(basename $(f)))/' >> $@${\n})
	tail -n +2 ontology/import.tsv | sed -e 's/$$/\timport/' >> $@

# Replace labels for any term with single quotes used in a C ROBOT template string
build/%-fixed.tsv: src/scripts/replace_labels.py build/index.tsv ontology/%.tsv
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
mro.xlsx: src/scripts/tsv2xlsx.py iedb/iedb.tsv $(source_files) iedb/iedb-manual.tsv ontology/rejected.tsv
	python3 $< $@ $(wordlist 2,100,$^)

update-tsv: update-tsv-files sort build/whitespace.tsv

# Update TSV files from Excel
.PHONY: update-tsv-files
update-tsv-files:
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
	python3 src/scripts/sort.py $(source_files)

# Sort IEDB table by IEDB ID
.PHONY: sort-iedb
sort-iedb: src/scripts/sort_iedb.py iedb/iedb.tsv
	python3 $^

# Check for whitespace during update-tsv step
.PRECIOUS: build/whitespace.tsv
build/whitespace.tsv: src/scripts/validation/detect_whitespace.py iedb/iedb.tsv iedb/iedb-manual.tsv $(source_files)
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
update-seqs: src/scripts/update_seqs.py ontology/chain.tsv build/hla.fasta build/mhc.fasta
	python3 $^

# refresh-seqs will overwrite existing seqs with new seqs
.PHONY: refresh-seqs
refresh-seqs: src/scripts/update_seqs.py ontology/chain.tsv build/hla.fasta build/mhc.fasta
	python3 $^ -o

# refresh-hla-seqs will overwrite existing HLA seqs with new seqs
.PHONY: refresh-hla-seqs
refresh-hla-seqs: src/scripts/update_seqs.py ontology/chain.tsv build/hla.fasta
	python3 $^ -o -H

build/hla_prot.fasta: | build
	curl -o $@ -L https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_prot.fasta

build/AlleleList.txt: | build
	curl -o $@ -L https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt

# TODO: update to not use index
.PHONY: update-mhcflurry-alleles
update-mhcflurry-alleles: src/scripts/alleles/update_alleles_mhcflurry.py ontology/chain-sequence.tsv ontology/chain.tsv ontology/molecule.tsv ontology/genetic-locus.tsv index.tsv build/mhc.fasta iedb/iedb.tsv
	python3 $^

### OWL Files

mro.owl: build/mro-import.owl build/index.tsv $(build_files) ontology/metadata.ttl | build/robot.jar
	$(ROBOT) template \
	--input $< \
	--merge-before \
	--template $(word 2,$^) \
	$(foreach i,$(build_files),--template $(i)) \
	reason \
	--reasoner ELK \
	--remove-redundant-subclass-axioms false \
	annotate \
	--ontology-iri "$(OBO)/mro.owl" \
	--version-iri "$(OBO)/mro/$(TODAY)/mro.owl" \
	--annotation owl:versionInfo "$(TODAY)" \
	--annotation-file ontology/metadata.ttl \
	--output $@

build/mro-import.owl: build/eco-import.ttl build/iao-import.ttl build/obi-import.ttl build/ro-import.ttl | build/robot.jar
	$(ROBOT) merge \
	$(foreach i,$^, --input $(i)) \
	annotate \
	--ontology-iri "$(OBO)/mro/$@" \
	extract \
	--method MIREOT \
	--upper-term "GO:0008150" \
	--upper-term "IAO:0000030" \
	--upper-term "OBI:1110128" \
	--upper-term "ECO:0000000" \
	--upper-term "BFO:0000040" \
	--lower-terms build/import.txt \
	remove \
	--term "CHEBI:23367" \
	--select "self descendants" \
	--exclude-term "PR:000000001" \
	--output $@

# fetch ontology dependencies
$(LIB)/%:
	mkdir -p $(LIB)
	cd $(LIB) && curl -LO "$(OBO)/$*"

UC = $(shell echo '$1' | tr '[:lower:]' '[:upper:]')

build/import.txt: ontology/import.tsv | build
	tail -n +2 $< | cut -f1 > $@

# Get all the OBI & GO terms plus RO:0000056, which is not in RO
build/obi.txt: build/import.txt
	sed '/^ECO/d' $< | sed '/^RO/d' | sed '/^IAO/d' > $@
	echo "RO:0000056" >> $@

# For each import, get just the terms in that namespace (OBI is a special case)
build/%.txt: build/import.txt
	sed -n '/$(call UC,$(notdir $(basename $@)))/p' $< > $@

# Create a RDFTab database for each import
# TODO: replace with LDTab
build/%.db: src/prefix.tsv $(LIB)/%.owl | build/rdftab
	rm -rf $@
	sqlite3 -cmd ".mode tabs" $@ ".import $< prefix"
	./build/rdftab $@ < $(word 2,$^)

PREDICATES := IAO:0000111 IAO:0000112 IAO:0000115 IAO:0000119 IAO:0000412 OBI:9991118 rdfs:label

# Extract with import module
# TODO: replace with gadget.extract
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
build/iedb.tsv: src/scripts/add_chain_accessions.py iedb/iedb.tsv ontology/molecule.tsv ontology/chain.tsv
	python3 $^ $@

# extended version for IEDB use
build/mro-iedb.owl: mro.owl build/iedb.tsv iedb/iedb-manual.tsv | build/robot.jar iedb
	$(ROBOT) template \
	--input $< \
	--template $(word 2,$^) \
	--template $(word 3,$^) \
	--merge-before \
	--output $@

build/.mro-tdb: build/mro-iedb.owl
	$(ROBOT) query --input $< \
	--create-tdb true \
	--tdb-directory $@

build/mhc_allele_restriction.csv: build/.mro-tdb src/queries/mhc_allele_restriction.rq | build/robot.jar
	$(ROBOT) query \
	--tdb-directory $< \
	--keep-tdb-mappings true \
	--select $(word 2,$^) $@

build/mhc_allele_restriction.tsv: src/scripts/clean.py build/mhc_allele_restriction.csv ontology/external.tsv | iedb
	python3 $^ > $@

build/molecule_export.tsv: src/scripts/export_molecule.py ontology/external.tsv ontology/molecule.tsv
	python3 $^ $@

build/ALLELE_FINDER_NAMES.csv: build/.mro-tdb src/queries/names.rq | build/robot.jar iedb
	$(ROBOT) query \
	--tdb-directory $< \
	--keep-tdb-mappings true \
	--select $(word 2,$^) $@.tmp \
	--format csv
	tail -n+2 $@.tmp | dos2unix > $@
	rm $@.tmp

build/ALLELE_FINDER_SEARCH.csv: build/.mro-tdb src/queries/search.rq | build/robot.jar iedb
	$(ROBOT) query \
	--tdb-directory $< \
	--keep-tdb-mappings true \
	--select $(word 2,$^) $@.tmp \
	--format csv
	tail -n+2 $@.tmp | dos2unix > $@
	rm $@.tmp

build/parents.csv: build/.mro-tdb src/queries/parents.rq | build/robot.jar
	$(ROBOT) query \
	--tdb-directory $< \
	--keep-tdb-mappings true \
	--select $(word 2,$^) $@

build/ALLELE_FINDER_TREE.csv: src/scripts/tree.py build/parents.csv | iedb
	python3 $^ --mode CSV > $@

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
# test: build/validation_errors_strict.tsv  TODO: new tests without index
test: build/report.csv
test: verify
test: build/mhc_allele_restriction_errors.tsv

# Python testing
.PHONY: pytest
pytest:
	py.test src/scripts/tree.py
	py.test src/scripts/synonyms.py

# Table diffs

DIFF_TABLES := build/diff/iedb.html $(foreach S,$(tables),build/diff/$(S).html) build/diff/iedb-manual.html

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

build/diff.html: src/scripts/diff.py src/scripts/diff.html $(DIFF_TABLES)
	python3 $< src/scripts/diff.html $(tables) > $@


### General

# Prepare products for testing & review
.PHONY: prepare
prepare: sort
prepare: clean
prepare: build/mro.db
prepare: update-seqs
prepare: update-iedb
prepare: build/diff.html

.PHONY: clean
clean:
	rm -rf mro.owl
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


### MRO Prototype

sql_tables = iedb iedb-manual $(tables) import rejected
sql_inputs = $(foreach t,$(sql_tables),build/tables/$(t).tsv)

# The tables for the SQL database do not include ROBOT template strings
.PHONY: copy_tables
copy_tables: $(sql_inputs)

build/tables/rejected.tsv: ontology/rejected.tsv | build/tables
	cp $< $@

build/tables/import.tsv: ontology/import.tsv | build/tables
	cp $< $@

build/tables/%.tsv: ontology/%.tsv | build/tables
	sed '2d' $< > $@

build/tables/iedb.tsv: iedb/iedb.tsv | build/tables
	sed '2d' $< > $@

build/tables/iedb-manual.tsv: iedb/iedb-manual.tsv | build/tables
	sed '2d' $< > $@

# Load all tables into SQLite database
build/mro-tables.db: src/scripts/load.py src/table.tsv src/column.tsv src/datatype.tsv $(sql_inputs)
	python3 src/scripts/load.py src/table.tsv $@

# Then add MRO using LDTab
.PHONY: load_mro
load_mro: mro.owl | build/ldtab.jar
	sqlite3 build/mro-tables.db "DROP TABLE IF EXISTS mro;"
	sqlite3 build/mro-tables.db "CREATE TABLE mro (assertion INT NOT NULL, retraction INT NOT NULL DEFAULT 0, graph TEXT NOT NULL, subject TEXT NOT NULL, predicate TEXT NOT NULL, object TEXT NOT NULL, datatype TEXT NOT NULL, annotation TEXT);"
	$(LDTAB) import --table mro build/mro-tables.db $<

# Delete existing database and reload (WARNING: this takes a long time!)
.PHONY: refresh_db
refresh_db:
	rm -rf build/mro-tables.db
	make load_mro

# Replace existing templates with data from database
.PHONY: refresh_templates
refresh_templates: $(foreach t,$(sql_tables),refresh_$(t))

refresh_iedb: src/scripts/table_2_template.py
	python3 $< build/mro-tables.db iedb iedb

refresh_iedb-manual: src/scripts/table_2_template.py
	python3 $< build/mro-tables.db iedb iedb-manual

refresh_%: src/scripts/table_2_template.py
	python3 $< build/mro-tables.db ontology $*
