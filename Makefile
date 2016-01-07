OBO = http://purl.obolibrary.org/obo
MRO = https://github.com/IEDB/MRO/raw/master/MRO_UNSTABLE.owl
LIB = lib

tables = external core loci haplotypes serotypes chains molecules haplotype-molecules serotype-molecules mutant-molecules sequences evidence-codes
TSVFiles = $(foreach o,$(tables),ontology/$(o).tsv)
templates = $(foreach i,$(TSVFiles),--template $(i))

# core
MRO_UNSTABLE.owl: mro-imports.owl index.tsv $(TSVFiles) ontology/annotations.ttl
	robot merge \
	--input mro-imports.owl \
	template \
	--prefix "MRO: $(MRO)#" \
	--prefix "REO: $(OBO)/REO_" \
	--template index.tsv \
	$(templates) \
	--merge-before \
	reason --reasoner HermiT \
	annotate \
	--ontology-iri "$(MRO)" \
	--version-iri "$(MRO)#$(shell date +%Y-%m-%d)" \
	--annotation-file ontology/annotations.ttl \
	--output $@

# extended version for IEDB use
MRO_UNSTABLE_IEDB.owl: MRO_UNSTABLE.owl ontology/iedb.tsv ontology/iedb-manual.tsv
	robot template \
	--prefix "MRO: $(MRO)#" \
	--input MRO_UNSTABLE.owl \
	--template ontology/iedb.tsv \
	--template ontology/iedb-manual.tsv \
	--merge-before \
	--output $@

# imports
mro-imports.owl: ontology/imports.txt $(LIB)/ro.owl $(LIB)/obi.owl $(LIB)/eco.owl
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

clean:
	rm -f mro-imports.owl MRO_UNSTABLE.owl MRO_UNSTABLE_IEDB.owl
