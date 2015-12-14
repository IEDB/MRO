# MHC Restriction Ontology (MRO)

<strong>WARNING: This is preliminary work in progress! Term identifiers are UNSTABLE, WILL NOT be maintained, and SHOULD NOT be used except for review!</strong>

[Download preliminary and unstable version of MRO.](https://github.com/IEDB/MRO/raw/master/MRO_UNSTABLE.owl)
We recommend [Protégé 5 beta](http://protege.stanford.edu/products.php#desktop-protege) for viewing MRO in OWL format.

MHC molecules form a highly diverse family of proteins that play a key role in cellular immune recognition. Over time, different techniques and terminologies have been developed to identify the specific type(s) of MHC molecule involved in a specific immune recognition context. No consistent nomenclature exists across different vertebrate species. To correctly represent MHC related data in [The Immune Epitope Database (IEDB)](http://www.iedb.org), we built upon a previously established MHC ontology ([MaHCO](http://www.bioinformatics.org/mahco/wiki/)) and created **MRO** to represent MHC molecules as they relate to immunological experiments. MRO models MHC protein chains from 16 species, deals with different approaches used to identify MHC, such as direct sequencing verses serotyping, relates engineered MHC molecules to naturally occurring ones, connects genetic loci, alleles, protein chains and multichain proteins, and establishes evidence codes for MHC restriction. Where available, this work is based on existing ontologies from the [OBO Foundry](http://obofoundry.org). We link to well-established MHC nomenclature resources of the international [ImMunoGeneTics information system (IMGT)](http://www.imgt.org) for human data and [The Immuno Polymorphism Database (IPD)](http://www.ebi.ac.uk/ipd) for non-human species.

This is **work in progress**. Here are some remaining tasks:

- model evidence codes
- have an expert review rat terms
- add sequences for the remaining chains, where possible
- review modelling with others ontology developers and domain experts
- add textual definitions and additional annotations

When the ontology is ready for wider use, we will migrate to stable term identifiers. The current term identifiers are UNSTABLE!


## Build Instructions

We use [ROBOT](https://github.com/ontodev/robot) and [GNU Make](https://www.gnu.org/software/make/) to build MRO from template files, applying the Quick Term Template (QTT) approach (see [Overcoming the ontology enrichment bottleneck with Quick Term Templates](http://dx.doi.org/10.3233/AO-2011-0086)).

We have one template file for the core (upper-level) terms:

- [mro-core.csv](mro-core.csv)

and then one file per branch of the ontology:

- [mro-loci.csv](mro-loci.csv)
- [mro-haplotypes.csv](mro-haplotypes.csv)
- [mro-serotypes.csv](mro-serotypes.csv)
- [mro-molecules.csv](mro-molecules.csv)
- [mro-chains.csv](mro-chains.csv)
- [mro-haplotype-molecules.csv](mro-haplotype-molecules.csv)
- [mro-serotype-molecules.csv](mro-serotype-molecules.csv)
- [mro-mutant-molecules.csv](mro-mutant-molecules.csv)
- [mro-sequences.csv](mro-sequences.csv)
- [mro-evidence-codes.csv](mro-evidence-codes.csv)

Two other files contain all remaining terms:

- [mro-obo.txt](mro-obo.txt) lists the external ontology terms that we import
- [mro-manual.owl](mro-manual.owl) contains object properties, annotation properties, and axioms that ROBOT cannot yet build from templates

With GNU Make and ROBOT installed, building is as simple as running:

	make

