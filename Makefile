OBO = http://purl.obolibrary.org/obo
MRO = https://github.com/IEDB/MRO/raw/master/MRO_UNSTABLE.owl
LIB = lib

generatedTables = loci haplotypes serotypes chains molecules haplotype-molecules serotype-molecules mutant-molecules
generatedCSVFiles = $(foreach o,$(generatedTables),mro-$(o).csv)
generatedTemplates = $(foreach i,$(generatedCSVFiles),--template $(i))

# core
MRO_UNSTABLE.owl: mro-manual.owl mro-obo.owl $(generatedCSVFiles) mro-sequences.csv mro-evidence-codes.csv
	robot merge \
	--input mro-manual.owl \
	--input mro-obo.owl \
	template \
	--prefix "MRO: $(MRO)#" \
	--template mro-core.csv \
	$(generatedTemplates) \
	--template mro-sequences.csv \
	--template mro-evidence-codes.csv \
	--merge-before \
	reason --reasoner HermiT \
	annotate \
	--version-iri "$(MRO)#$(shell date +%Y-%m-%d)" \
	--output MRO_UNSTABLE.owl

# extended version for IEDB use
MRO_UNSTABLE_IEDB.owl: MRO_UNSTABLE.owl mro-iedb.csv mro-iedb-manual.csv
	robot template \
	--prefix "MRO: $(MRO)#" \
	--input MRO_UNSTABLE.owl \
	--template mro-iedb.csv \
	--template mro-iedb-manual.csv \
	--merge-before \
	--output MRO_UNSTABLE_IEDB.owl

# imports
mro-obo.owl: mro-obo.txt $(LIB)/ro.owl $(LIB)/obi.owl $(LIB)/eco.owl
	robot merge \
	--input $(LIB)/ro.owl \
	--input $(LIB)/obi.owl \
	--input $(LIB)/eco.owl \
	extract \
	--upper-term "GO:0008150" \
	--upper-term "IAO:0000030" \
	--upper-term "OBI:1110128" \
	--upper-term "ECO:0000000" \
	--upper-term "BFO:0000040" \
	--upper-term "PR:000000001" \
	--lower-terms mro-obo.txt \
	--output mro-obo.owl

# fetch ontology dependencies
$(LIB)/%:
	mkdir -p $(LIB)
	cd $(LIB) && curl -LO "$(OBO)/$*"

clean:
	rm -f mro-obo.owl MRO_UNSTABLE.owl MRO_UNSTABLE_IEDB.owl
