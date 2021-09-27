import csv
import os
import sys
from subprocess import Popen, PIPE
import json

class Species:
    """Species object contains metadata required to have flexible logic regarding MHC
    nomenclature. This includes things such as MHC Code (Mafa, Mamu, etc.), taxon,
    class I alleles, class II alpha chains and their corresponding beta chains.

    This is all loaded through locus_data.json, which can be extended to any new
    species or any changes in nomenclature logic (Mafa-AG -> Mafa-AG1, AG2)"""

    def __init__(self, data):
        self.mhc_code = data["mhc_code"]
        self.taxon = data["taxon"]
        self.class_I = data["class_I"]
        self.pairing = data["pairing"]
        self.pairing_rev = data["pairing_rev"]

    def get_class_I(self):
        """Return a set containing strictly class I loci"""
        return set(self.class_I)

    def get_class_II(self):
        """Return a set containing strictly class II loci"""
        class_II = set()
        # This is a dict with key: value (allele 1: allele 2) for alpha/beta chain pair
        for allele_key in self.pairing.keys():
            # Allele key is first allele
            class_II.add(allele_key)
            # Value from pairing is second allele
            class_II.add(self.pairing[allele_key])


def get_MRO_seqs():
    """Get both the MRO alleles and the MRO generic chains"""
    mro_alleles = set()
    mro_gen_chains = set()
    with open(sys.argv[2]) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            allele = row["Label"].split(" ")[0]
            allele = (":").join(allele.split(":")[:2])
            mro_alleles.add(allele)
            # If we don't have a * in an allele it is a generic chain i.e. "Mafa-A1 Chain"
            if "*" not in row["Label"]:
                mro_gen_chains.add(row["Label"])

    return mro_alleles, mro_gen_chains


def get_IPD_seqs():
    """Gets sequences from IPD-MHC database"""
    # Create a dictionary mapping the IMGT accession to protein sequence
    seqs = {}
    allele_names = {}
    with open(sys.argv[6]) as fasta:
        accession = None
        seq = ""
        allele = ""
        for line in fasta:
            if line.startswith(">"):
                if accession:
                    seqs[accession] = seq
                    allele_names[allele] = accession

                accession = None
                seq = ""
                allele = ""

                # Match the accession
                if line.startswith(">IPD-MHC"):
                    accession = line.split(" ")[0][9:]
                    allele = line.split(" ")[1]
                    allele = (":").join(allele.split(":")[:4])
            else:
                seq += line.strip()
        seqs[accession] = seq
        allele_names[allele] = accession
    return seqs, allele_names


def get_species_metadata():
    """Reads in a JSON file with metadata about the various species in MRO
    This JSON file is easily extensible to any new species or new loci
    """
    # Load the JSON data
    with open(sys.argv[8], "r") as fh:
        data = json.load(fh)

    class_I = set()
    class_II = set()

    species_map = {}
    # Generate the map of mhc_code: Species objects
    for species_data in data["species"]:
        species_map[species_data["mhc_code"]] = Species(species_data)

        # Get a set of class I and class II alleles for protein complex logic
        for allele in species_data["class_I"]:
            class_I.add(f"{species_data['mhc_code']}-{allele}")
        for alpha_allele in species_data["pairing"].keys():
            # Get the alpha chain allele which maps to the beta chain allele
            class_II.add(f"{species_data['mhc_code']}-{alpha_allele}")
            # species_data['pairing'][alpha_allele] just gives the mapped beta chain allele
            class_II.add(
                f"{species_data['mhc_code']}-{species_data['pairing'][alpha_allele]}"
            )

    return species_map, class_I, class_II


def get_missing_alleles(ipd_alleles, mro_alleles, classI, classII):
    """Get missing alleles by comparing IPD-MHC alleles and what is in MRO

    Args:
        ipd_alleles: key value pair of allele: sequence
        mro_alleles: list of mro alleles
    """
    # Load the args into sets for easy comparisons
    ipd_alleles_set = set(ipd_alleles.keys())
    mro_alleles_set = set(mro_alleles)

    # Find items in ipd-mhc not in MRO
    missing_alleles = ipd_alleles_set.difference(mro_alleles_set)
    # Make sure we aren't processing any HLA alleles by accident
    missing_alleles = {x for x in missing_alleles if "HLA" not in x}
    missing_alleles = {x for x in missing_alleles if x.split("*")[0] in classI or x.split("*")[0] in classII}

    return missing_alleles


def get_MRO_loci():
    """Get the current loci from genetic-locus.tsv"""
    mro_loci = set()
    with open(sys.argv[4]) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            locus = row["Label"]
            mro_loci.add(locus.split(" ")[0])

    return mro_loci


def create_new_loci(species_map, curr_MRO_loci):
    """Update genetic-locus.tsv with missing loci"""
    new_loci_rows = set()
    missing_loci = set()

    # Species map contains mhc_code: Species object
    for species in species_map.values():

        # Class I loci
        for classI in species.get_class_I():
            allele = species.mhc_code + "-" + classI
            if allele not in curr_MRO_loci:
                missing_loci.add(allele)
                mhc_class = "MHC class I"
                label = f"{allele} locus"
                synonyms = ""
                class_type = "subclass"
                parent = f"{species.taxon} {mhc_class} locus"
                new_loci_rows.add((label, synonyms, class_type, parent, ""))

        # Class II loci
        for classII in species.pairing.keys():
            allele1 = species.mhc_code + "-" + classII
            allele2 = species.mhc_code + "-" + species.pairing[classII]

            # Build the tuple for the first allele
            if allele1 not in curr_MRO_loci:
                missing_loci.add(allele1)
                mhc_class = "MHC class II"
                label = f"{allele1} locus"
                synonyms = ""
                class_type = "subclass"
                parent = f"{species.taxon} {mhc_class} locus"
                new_loci_rows.add((label, synonyms, class_type, parent, ""))

            # Build the tuple for the second allele
            if allele2 not in curr_MRO_loci:
                missing_loci.add(allele2)
                mhc_class = "MHC class II"
                label = f"{allele2} locus"
                synonyms = ""
                class_type = "subclass"
                parent = f"{species.taxon} {mhc_class} locus"
                new_loci_rows.add((label, synonyms, class_type, parent, ""))

    with open(sys.argv[4], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in new_loci_rows:
            writer.writerow(tup)

    return missing_loci


def create_new_chains(new_alleles, mro_gen_chains, ipd_seqs, allele_map):
    # Declare sets to make sure there are no duplicate terms
    missing_gen_chains = set()
    missing_gen_chain_rows = set()

    missing_chain_rows = set()
    missing_chain_seq_rows = set()

    for allele in new_alleles:
        gene = allele.split("*")[0]
        # Check if we have generic chain i.e. "Mafa-A1 chain" in template
        if gene + " chain" not in mro_gen_chains:
            missing_gen_chain_rows.add(
                (gene + " chain", "", "equivalent", "protein", gene + " locus")
            )
            missing_gen_chains.add(gene)
        # Now move on to specific allele chain i.e. "Mafa-A1*01:01 chain"
        try:
            # Add to the chain.tsv template
            accession = allele_map[allele]
            source = "IPD"
            label = f"{allele} chain"
            resource_name = allele
            seq = ipd_seqs[accession]
            missing_chain_seq_rows.add((label, resource_name, source, accession, seq))

            # Add to the chain-sequence.tsv template
            synonyms = ""
            class_type = "subclass"
            parent = f"{gene} chain"
            missing_chain_rows.add((label, synonyms, class_type, parent, "", ""))

        except Exception:
            print(f"Allele {allele} has no sequence in IPD")

    # Write out new rows to chain.tsv template
    with open(sys.argv[2], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in missing_chain_rows:
            writer.writerow(tup)
        for tup in missing_gen_chain_rows:
            writer.writerow(tup)

    # Write out new rows to chain-sequence.tsv template
    with open(sys.argv[1], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in missing_chain_seq_rows:
            writer.writerow(tup)

    return missing_gen_chains


def create_generic_protein(mro_proteins, gene, species, classI):
    if gene not in mro_proteins:
        label = f"{gene} protein complex"
        iedb_name = f"{gene}"
        synonym = ""
        if classI:
            class_type = "equivalent"
            parent = "MHC class I protein complex"
            alpha_chain = f"{gene} chain"
            beta_chain = f"Beta-2-microglobulin"
            restrict_lvl = "locus"
        else:
            loci = gene.split("-")[1]
            parent = "MHC class II protein complex"
            alpha_chain = f"{gene} chain"
            try:
                beta_chain = f"{species.pairing[loci]} chain"
            except:
                beta_chain = f"{species.pairing_rev[loci]} chain"
            restrict_lvl = "partial molecule"

        new_generic_protein = (
            label,
            iedb_name,
            synonym,
            restrict_lvl,
            class_type,
            parent,
            species.taxon,
            alpha_chain,
            beta_chain,
            "",
            "",
        )
        mro_proteins.add(f"{gene} chain")

        with open(sys.argv[3], "a+") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(new_generic_protein)
        
        return label
    else:
        return None


def create_protein(missing_proteins, classI):
    new_molecules = set()
    new_molecule_labels = set()
    for tup in missing_proteins:
        allele = tup[0]
        species = tup[1]
        label = f"{allele} protein complex"
        iedb_name = f"{allele}"
        synonym = ""
        class_type = "equivalent"
        if classI:
            restrict_lvl = "complete molecule"
            parent = "MHC class I protein complex"
            alpha_chain = f"{allele} chain"
            beta_chain = "Beta-2-microglobulin"
        else:
            loci = allele.split("-")[1].split("*")[0]
            restrict_lvl = "partial molecule"
            parent = "MHC class II protein complex"
            alpha_chain = f"{allele} chain"
            try:
                beta_chain = f"{species.mhc_code}-{species.pairing[loci]} chain"
            except:
                beta_chain = f"{species.mhc_code}-{species.pairing_rev[loci]} chain"
        new_protein = (
                label,
                iedb_name,
                synonym,
                restrict_lvl,
                class_type,
                parent,
                species.taxon,
                alpha_chain,
                beta_chain,
                "",
                ""
        )
        new_molecules.add(new_protein)
        new_molecule_labels.add(f"{allele} protein complex")
    
    with open(sys.argv[3], "a+") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for new_molecule in new_molecules:
            writer.writerow(new_molecule)

    return new_molecule_labels


def create_new_molecules(missing_alleles, class_I, class_II, species_map):
    mro_proteins = set()
    new_molecules = set()
    new_generic_molecules = set()
    with open(sys.argv[3]) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            mro_protein = row["IEDB Label"]
            if "/" in mro_protein:
                continue
            else:
                mro_proteins.add(mro_protein)

    new_class_I = set()
    new_class_II = set()
    for allele in missing_alleles:
        try:
            mhc_code = allele.split("-")[0]
            gene = allele.split("*")[0]
            species = species_map[mhc_code]
            if gene in class_I:
                gen_prot = create_generic_protein(mro_proteins, gene, species, True)
                new_class_I.add((allele, species))
                if gen_prot != None:
                    new_generic_molecules.add(gen_prot)
            if gene in class_II:
                gen_prot = create_generic_protein(mro_proteins, gene, species, False)
                new_class_II.add((allele, species))
                if gen_prot != None:
                    new_generic_molecules.add(gen_prot)
        except:
            print(f"{gene} not available in locus-data.json, please add it.")

    for x in list(create_protein(new_class_I, True)):
        if x not in mro_proteins:
            new_molecules.add(x)
    for y in list(create_protein(new_class_II, False)):
        if y not in mro_proteins:
            new_molecules.add(y)

    return new_molecules, new_generic_molecules

def update_index(missing_alleles, missing_molecules, missing_loci, new_gen_chains, new_gen_molecules):
    mro_labels = set()
    with open(sys.argv[5], "r") as fh:
        for _ in range(2):
            next(fh)
        for line in fh:
            curr_mro_id = line.rstrip().split("\t")[0]
            curr_mro_label = line.rstrip().split("\t")[1]
            mro_labels.add(curr_mro_label)
            curr_mro_num = int(curr_mro_id.split(":")[1])

    # Next few blocks will iterate MRO ID and pad left with 0s to 7 digits
    new_tups = []
    for chain in new_gen_chains:
        chain_name = f"{chain} chain"
        if chain_name not in mro_labels:
            mro_id = f"MRO:{str(curr_mro_num + 1).zfill(7)}"
            curr_mro_num += 1
            new_tups.append((mro_id, chain_name, "owl:Class"))

    for locus in missing_loci:
        locus_name = f"{locus} locus"
        if locus_name not in mro_labels:
            mro_id = f"MRO:{str(curr_mro_num + 1).zfill(7)}"
            curr_mro_num += 1
            new_tups.append((mro_id, locus_name, "owl:Class"))

    for allele in missing_alleles:
        chain_name = f"{allele} chain"
        if chain_name not in mro_labels:
            mro_id = f"MRO:{str(curr_mro_num + 1).zfill(7)}"
            curr_mro_num += 1
            new_tups.append((mro_id, chain_name, "owl:Class"))

    for molecule in missing_molecules:
        molecule_name = f"{molecule}"
        if molecule_name not in mro_labels:
            mro_id = f"MRO:{str(curr_mro_num + 1).zfill(7)}"
            curr_mro_num += 1
            new_tups.append((mro_id, molecule_name, "owl:Class"))

    for molecule in new_gen_molecules:
        molecule_name = f"{molecule}"
        if molecule_name not in mro_labels:
            mro_id = f"MRO:{str(curr_mro_num + 1).zfill(7)}"
            curr_mro_num += 1
            new_tups.append((mro_id, molecule_name, "owl:Class"))

    with open(sys.argv[5], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in new_tups:
            writer.writerow(tup)


def update_IEDB_tab(missing_molecules):
    iedb_labels = set()
    with open(sys.argv[7], "r") as fh:
        for _ in range(2):
            next(fh)
        for line in fh:
            curr_iedb_label = line.rstrip().split("\t")[0]
            iedb_labels.add(curr_iedb_label)
    
    # Get current IEDB id
    with open(sys.argv[7]) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
        final_row = rows[-1]
        curr_iedb_id = int(final_row["IEDB ID"])

    # Write out results
    with open(sys.argv[7], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for molecule in missing_molecules:
            if molecule not in iedb_labels:
                locus = molecule.split("-")[1].split("*")[0]
                curr_iedb_id += 1
                if "DQ" in locus:
                    tup = (molecule, curr_iedb_id, "DQ", "", "")
                elif "DRB" in locus:
                    tup = (molecule, curr_iedb_id, "DRB", "", "")
                else:
                    tup = (molecule, curr_iedb_id, locus, "", "")
                writer.writerow(tup)


def update_MRO():
    """Primary driver for updating MRO with IPD-MHC sequences

    Works down a hierarchy starting with missing alleles, then
    creating the templates for new loci, chains, and molecules

    Finish by adding to index and IEDB template"""

    species_map, class_I, class_II = get_species_metadata()

    # Get new alleles
    mro_alleles, mro_gen_chains = get_MRO_seqs()
    ipd_alleles, ipd_allele_map = get_IPD_seqs()
    
    new_alleles = get_missing_alleles(ipd_allele_map, mro_alleles, class_I, class_II)

    # Get new and current loci
    curr_MRO_loci = get_MRO_loci()

    new_loci = create_new_loci(species_map, curr_MRO_loci)

    # create and return new chains (gen + specific)
    new_gen_chains = create_new_chains(
        new_alleles, mro_gen_chains, ipd_alleles, ipd_allele_map
    )

    # create and return new molecules (gen + specific)
    new_molecules, new_gen_molecules = create_new_molecules(new_alleles, class_I, class_II, species_map)

    # Update index
    update_index(new_alleles, new_molecules, new_loci, new_gen_chains, new_gen_molecules)

    # Update IEDB
    update_IEDB_tab(new_molecules)

update_MRO()