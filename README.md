# MHC Restriction Ontology (MRO)

[![Build Status](https://travis-ci.com/IEDB/MRO.svg?branch=master)](https://travis-ci.com/IEDB/MRO)

[**Download the latest release of MRO.**](http://purl.obolibrary.org/obo/mro.owl)

[Download the development version of MRO.](mro.owl?raw=true)

We recommend [Protégé 5 beta](http://protege.stanford.edu/products.php#desktop-protege) for viewing MRO in OWL format.

MHC molecules form a highly diverse family of proteins that play a key role in cellular immune recognition. Over time, different techniques and terminologies have been developed to identify the specific type(s) of MHC molecule involved in a specific immune recognition context. No consistent nomenclature exists across different vertebrate species. To correctly represent MHC related data in [The Immune Epitope Database (IEDB)](http://www.iedb.org), we built upon a previously established MHC ontology ([MaHCO](http://www.bioinformatics.org/mahco/wiki/)) and created **MRO** to represent MHC molecules as they relate to immunological experiments. MRO models MHC protein chains from 16 species, deals with different approaches used to identify MHC, such as direct sequencing verses serotyping, relates engineered MHC molecules to naturally occurring ones, connects genetic loci, alleles, protein chains and multichain proteins, and establishes evidence codes for MHC restriction. Where available, this work is based on existing ontologies from the [OBO Foundry](http://obofoundry.org). We link to well-established MHC nomenclature resources of the international [ImMunoGeneTics information system (IMGT)](http://www.imgt.org) for human data and [The Immuno Polymorphism Database (IPD)](http://www.ebi.ac.uk/ipd) for non-human species.

Vita et al., [An Ontology for Major Histocompatibility Restriction](http://www.jbiomedsem.com/content/7/1/1), Journal of Biomedical Semantics, 2016.

This is **work in progress**. Here are some remaining tasks:

- model evidence codes
- have an expert review rat terms
- add sequences for the remaining chains, where possible
- review modelling with others ontology developers and domain experts
- add textual definitions and additional annotations


## Build Instructions

We use [ROBOT](https://github.com/ontodev/robot) and [GNU Make](https://www.gnu.org/software/make/) to build MRO from template files, applying the Quick Term Template (QTT) approach (see [Overcoming the ontology enrichment bottleneck with Quick Term Templates](http://dx.doi.org/10.3233/AO-2011-0086)).

The [index.tsv](index.tsv) contains the master list of MRO terms. The [ontology/](ontology/) directory contains template files, with one file for the core (upper-level) terms:

- [core.tsv](ontology/core.tsv)

and then one file per branch of the ontology:

- [genetic-locus.tsv](ontology/genetic-locus.tsv)
- [haplotype.tsv](ontology/haplotype.tsv)
- [serotype.tsv](ontology/serotype.tsv)
- [chain.tsv](ontology/chain.tsv)
- [molecule.tsv](ontology/molecule.tsv)
- [haplotype-molecule.tsv](ontology/haplotype-molecule.tsv)
- [serotype-molecule.tsv](ontology/serotype-molecule.tsv)
- [mutant-molecule.tsv](ontology/mutant-molecule.tsv)
- [evidence.tsv](ontology/evidence.tsv)
- [chain-sequence.tsv](ontology/chain-sequence.tsv)

Two other files contain all remaining terms:

- [import.txt](ontology/import.txt) lists the external ontology terms that we import
- [external.tsv](ontology/external.tsv) declares some external ontology terms

With GNU Make and ROBOT installed, building is as simple as running:

	make

