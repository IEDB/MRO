OBO = http://purl.obolibrary.org/obo
LIB = lib

tables = external core genetic-locus haplotype serotype chain molecule haplotype-molecule serotype-molecule mutant-molecule evidence chain-sequence
source_files = $(foreach o,$(tables),ontology/$(o).tsv)
build_files = $(foreach o,$(tables),build/$(o).tsv)
templates = $(foreach i,$(build_files),--template $(i))

build:
	mkdir -p $@

build/%.tsv: ontology/%.tsv | build
	cp $< $@
build/molecule.tsv: src/synonyms.py ontology/molecule.tsv | build
	$^ > $@
build/haplotype-molecule.tsv: src/synonyms.py ontology/haplotype-molecule.tsv | build
	$^ > $@
build/serotype-molecule.tsv: src/synonyms.py ontology/serotype-molecule.tsv | build
	$^ > $@
build/mutant-molecule.tsv: src/synonyms.py ontology/mutant-molecule.tsv | build
	$^ > $@


# Represent tables in Excel

mro.xlsx: src/tsv2xlsx.py index.tsv iedb/iedb.tsv ontology/genetic-locus.tsv ontology/haplotype.tsv ontology/serotype.tsv ontology/chain.tsv ontology/chain-sequence.tsv ontology/molecule.tsv ontology/haplotype-molecule.tsv ontology/serotype-molecule.tsv ontology/mutant-molecule.tsv ontology/core.tsv ontology/external.tsv iedb/iedb-manual.tsv ontology/evidence.tsv
	$< $@ $(wordlist 2,100,$^)

.PHONY: update-tsv
update-tsv:
	src/xlsx2tsv.py mro.xlsx index > index.tsv
	src/xlsx2tsv.py mro.xlsx iedb > iedb/iedb.tsv
	src/xlsx2tsv.py mro.xlsx iedb-manual > iedb/iedb-manual.tsv
	$(foreach t,$(tables),src/xlsx2tsv.py mro.xlsx $(t) > ontology/$(t).tsv;)
	src/sort.py $(source_files)



# core
mro.owl: mro-import.owl index.tsv $(build_files) ontology/metadata.ttl
	robot merge \
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

# import
mro-import.owl: ontology/import.txt $(LIB)/ro.owl $(LIB)/obi.owl $(LIB)/eco.owl
	robot merge \
	--input $(LIB)/eco.owl \
	--input $(LIB)/obi.owl \
	--input $(LIB)/ro.owl \
	extract \
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


# Utilities

.PHONY: sort
sort:
	src/sort.py $(source_files)


# generate files for IEDB

iedb:
	mkdir -p $@

# extended version for IEDB use
iedb/mro-iedb.owl: mro.owl iedb/iedb.tsv iedb/iedb-manual.tsv | iedb
	robot template \
	--prefix "MRO: $(OBO)/MRO_" \
	--input $< \
	--template $(word 2,$^) \
	--template $(word 3,$^) \
	--merge-before \
	--output $@

build/mhc_allele_restriction.csv: iedb/mro-iedb.owl src/mhc_allele_restriction.rq | build
	robot query --input $(word 1,$^) --select $(word 2,$^) $@

iedb/mhc_allele_restriction.tsv: src/clean.py build/mhc_allele_restriction.csv | iedb
	python3 $^ \
	> $@

iedb/ALLELE_FINDER_NAMES.csv: iedb/mro-iedb.owl src/names.rq | iedb
	robot query --input $(word 1,$^) --select $(word 2,$^) $@.tmp
	tail -n+2 $@.tmp | dos2unix > $@
	rm $@.tmp

iedb/ALLELE_FINDER_SEARCH.csv: iedb/mro-iedb.owl src/search.rq | iedb
	robot query --input $(word 1,$^) --select $(word 2,$^) $@.tmp
	tail -n+2 $@.tmp | dos2unix > $@
	rm $@.tmp

build/parents.csv: iedb/mro-iedb.owl src/parents.rq | build
	robot query --input $(word 1,$^) --select $(word 2,$^) $@

iedb/ALLELE_FINDER_TREE.csv: src/tree.py build/parents.csv | iedb
	$^ --mode CSV > $@

build/tree.json: src/tree.py build/parents.csv | build
	$^ --mode JSON > $@

build/full_tree.json: src/tree.py build/full_tree.csv | build
	$^ --mode JSON > $@


IEDB_TARGETS := iedb/mro-iedb.owl iedb/mhc_allele_restriction.tsv iedb/ALLELE_FINDER_NAMES.csv iedb/ALLELE_FINDER_SEARCH.csv iedb/ALLELE_FINDER_TREE.csv
.PHONY: update-iedb
update-iedb: $(IEDB_TARGETS)


.PHONY: test
test:
	py.test src/tree.py
	py.test src/synonyms.py

.PHONY: clean
clean:
	rm -rf build
	rm -f mro.owl mro-import.owl
	rm -f $(IEDB_TARGETS)
