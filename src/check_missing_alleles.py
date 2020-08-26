import csv
import os
from subprocess import Popen, PIPE


def get_IMGT_data():
    """Uses wget to retrieve the latest IMGT alleles + sequences from their github
    These get deleted later after operations are complete
    """
    imgt_prot_URL = (
        "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_prot.fasta"
    )
    imgt_allele_URL = (
        "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt"
    )
    process = Popen(["wget", imgt_prot_URL, imgt_allele_URL], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()


def update_chain_sequence():
    """Updates the chain-sequence.tsv with missing alleles from IMGT
    
    Returns:
        A list of missing alleles ['A*02:01', 'B*02:01']
    """

    mro_alleles = set()
    with open("../ontology/chain-sequence.tsv") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            if "HLA" in row["Label"]:
                mro_alleles.add(row["Label"].split(" ")[0][4:])

    # Create a dictionary mapping the IMGT accession to protein sequence
    seqs = {}
    with open("./hla_prot.fasta") as fasta:
        accession = None
        seq = ""
        for line in fasta:
            if line.startswith(">"):
                if accession:
                    seqs[accession] = seq
                accession = None
                seq = ""

                # Match the accession
                if line.startswith(">HLA:HLA"):
                    accession = line[5:13]
                else:
                    print("Bad accession:", line)
            else:
                seq += line.strip()
        seqs[accession] = seq

    imgt_dict = {}

    # Get the mapping between allele name and IMGT accession
    with open("./Allelelist.txt") as fh:
        # Skip IMGT header
        for _ in range(6):
            next(fh)
        rows = csv.DictReader(fh, delimiter=",")
        for row in rows:
            allele = row["Allele"]
            allele = (":").join(allele.split(":")[:2])
            if allele not in imgt_dict and "N" not in allele and "Q" not in allele:
                imgt_dict[allele] = row["AlleleID"]

    # Find differences between IMGT and MRO alleles
    imgt_alleles = set(imgt_dict.keys())
    missing_alleles = imgt_alleles.difference(mro_alleles)

    missing_allele_rows = set()
    for allele in missing_alleles:
        try:
            accession = imgt_dict[allele]
            source = "IMGT/HLA"
            label = f"HLA-{allele} chain"
            resource_name = allele
            seq = seqs[imgt_dict[allele]]
            missing_allele_rows.add((label, resource_name, source, accession, seq))
        except:
            continue

    with open("../ontology/chain-sequence.tsv", "a+") as outfile:
        for tup in missing_allele_rows:
            outfile.write(("\t").join(tup) + "\n")

    return missing_alleles


def update_chain(missing_alleles):
    """Update the chain.tsv table to add new chains from IMGT

    Returns:
        A list of missing genes ['DRB', 'Y']
    """

    missing_genes = set()
    for allele in missing_alleles:
        gene = allele.split("*")[0]
        missing_genes.add(gene)

    mro_chains = set()
    mro_genes = set()
    with open("../ontology/chain.tsv") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            label = row["Label"]
            if "HLA" in label:
                if "*" not in label:
                    gene = label.split(" ")[0][4:]
                    mro_genes.add(gene)
                else:
                    allele = label.split(" ")[0][4:]
                    mro_chains.add(allele)

    new_genes = missing_genes.difference(mro_genes)
    new_gene_tups = set()
    for gene in new_genes:
        label = f"HLA-{gene} chain"
        synonyms = ""
        class_type = "equivalent"
        parent = "protein"
        gene = f"HLA-{gene} locus"
        new_gene_tups.add((label, synonyms, class_type, parent, gene))

    new_chains = missing_alleles.difference(mro_chains)
    new_chain_tups = set()
    for allele in new_chains:
        label = f"HLA-{allele} chain"
        synonyms = ""
        class_type = "subclass"
        # This removes the specific allele # (generic chain)
        parent = f"HLA-{allele.split('*')[0]} chain"
        gene = ""
        new_chain_tups.add((label, synonyms, class_type, parent, gene))

    with open("../ontology/chain.tsv", "a+") as fh:
        for tup in new_gene_tups:
            fh.write(("\t").join(tup) + "\n")
        for tup in new_chain_tups:
            fh.write(("\t").join(tup) + "\n")

    return new_genes


def update_locus(missing_genes):
    """Update genetic-locus.tsv with any missing genes
    
    Returns:
        A list of missing genetic loci ['DRB', 'Y']
    """
    class_I = ["A", "B", "C"]
    class_II = ["DP", "DM", "DOA", "DOB", "DQ", "DR"]

    mro_loci = set()
    with open("../ontology/genetic-locus.tsv") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            locus = row["Label"]
            if "HLA" in locus:
                mro_loci.add(locus.split(" ")[0][4:])

    new_loci_rows = set()
    new_loci = missing_genes.difference(mro_loci)
    for locus in new_loci:
        # Fast and dirty way to check if gene is classI, classII
        mhc_class = ""
        for prefix in class_I:
            if prefix in locus:
                mhc_class = "MHC class I"
        for prefix in class_II:
            if prefix in locus:
                mhc_class = "MHC class II"
        if mhc_class == "":
            mhc_class = "non-classical MHC"

        label = f"HLA-{locus} locus"
        synonyms = ""
        class_type = "subclass"
        parent = f"human {mhc_class} locus"
        new_loci_rows.add((label, synonyms, class_type, parent))

    with open("../ontology/genetic-locus.tsv", "a+") as fh:
        for row in new_loci_rows:
            fh.write(("\t").join(row) + "\n")

    return new_loci


def update_index(missing_alleles, missing_genes, missing_loci, missing_molecules):
    """Adds new entries from IMGT to index.tsv to allow ROBOT to build owl file
    """
    mro_labels = set()
    with open("../index.tsv") as fh:
        for _ in range(2):
            next(fh)
        for line in fh:
            curr_mro_id = line.rstrip().split("\t")[0]
            curr_mro_label = line.rstrip().split("\t")[1]
            mro_labels.add(curr_mro_label)
            curr_mro_num = int(curr_mro_id.split(":")[1])

    # Next few blocks will iterate MRO ID and pad left with 0s to 7 digits
    new_tups = []
    for allele in missing_alleles:
        chain_name = f"HLA-{allele} chain"
        if chain_name not in mro_labels:
            mro_id = f"MRO:{str(curr_mro_num + 1).zfill(7)}"
            curr_mro_num += 1
            new_tups.append((mro_id, chain_name, "owl:Class"))

    for gene in missing_genes:
        gene_name = f"HLA-{gene} chain"
        if gene_name not in mro_labels:
            mro_id = f"MRO:{str(curr_mro_num + 1).zfill(7)}"
            curr_mro_num += 1
            new_tups.append((mro_id, gene_name, "owl:Class"))

    for locus in missing_loci:
        locus_name = f"HLA-{locus} locus"
        if locus_name not in mro_labels:
            mro_id = f"MRO:{str(curr_mro_num + 1).zfill(7)}"
            curr_mro_num += 1
            new_tups.append((mro_id, locus_name, "owl:Class"))

    for molecule in missing_molecules:
        molecule_name = f"HLA-{molecule} protein complex"
        if molecule_name not in mro_labels:
            mro_id = f"MRO:{str(curr_mro_num + 1).zfill(7)}"
            curr_mro_num += 1
            new_tups.append((mro_id, molecule_name, "owl:Class"))

    with open("../index.tsv", "a+") as fh:
        for tup in new_tups:
            fh.write(("\t").join(tup) + "\n")

def create_classI_prot(missing_chainI):
    new_classI_molecules = set()
    with open("../ontology/molecule.tsv", "a+") as fh:
        for allele in missing_chainI:
            label = f"HLA-{allele} protein complex"
            iedb_name = f"HLA-{allele}"
            synonym = ""
            restrict_lvl = "complete molecule"
            class_type = "equivalent"
            parent = "MHC class I protein complex"
            taxon = "human"
            alpha_chain = f"HLA-{allele} chain"
            beta_chain = "Beta-2-microglobulin"
            fh.write(("\t").join([label, iedb_name, synonym, restrict_lvl, class_type, parent, taxon, alpha_chain, beta_chain]) + "\n")
            new_classI_molecules.add(f"HLA-{allele} protein complex")
    return new_classI_molecules

def create_classII_pairing(allele):
    """Simple helper class to create a pairing of alpha and beta chain"""
    pairing = {"DQA1": "DQB1", "DRA": ["DRB1", "DRB3", "DRB4", "DRB5"], "DPA1": "DPB1"}
    pairing_rev = {"DQB1": "DQA1", "DRB1": "DRA", "DRB3": "DRA", "DRB4": "DRA", "DRB5": "DRA", "DPB1": "DPA1"}
    gene = allele.split("*")[0]
    if gene in pairing:
        pair = f"{allele}/{pairing[gene]}"
    else:
        pair = f"{pairing_rev[gene]}/{allele}"
    
    return pair

def create_classII_prot(missing_chainII):
    new_classII_molecules = set()
    with open("../ontology/molecule.tsv", "a+") as fh:
        for allele in missing_chainII:
            pair = create_classII_pairing(allele)
            alpha_gene = pair.split("/")[0]
            beta_gene = pair.split("/")[1]
            
            # First create a complete molecule entry
            label = f"HLA-{pair} protein complex"
            iedb_name = f"HLA-{pair}"
            synonym = ""
            restrict_lvl = "complete molecule"
            class_type = "equivalent"
            parent = "MHC class II protein complex"
            taxon = "human"
            alpha_chain = f"HLA-{alpha_gene} chain"
            beta_chain = f"HLA-{beta_gene} chain"
            
            fh.write(("\t").join([label, iedb_name, synonym, restrict_lvl, class_type, parent, taxon, alpha_chain, beta_chain]) + "\n")
            new_classII_molecules.add(f"HLA-{pair} protein complex")

            # Now partial molecule entry
            label = f"HLA-{allele} protein complex"
            iedb_name = f"HLA-{allele}"
            synonym = ""
            restrict_lvl = "partial molecule"
            class_type = "equivalent"
            parent = "MHC class II protein complex"
            taxon = "human"
            alpha_chain = f"HLA-{alpha_gene} chain"
            beta_chain = f"HLA-{beta_gene} chain"
            
            fh.write(("\t").join([label, iedb_name, synonym, restrict_lvl, class_type, parent, taxon, alpha_chain, beta_chain]) + "\n")
            new_classII_molecules.add(f"HLA-{allele} protein complex")
    
    return new_classII_molecules

def update_molecules(missing_alleles):
    class_I_genes = ["A", "B", "C"]
    class_II_genes = ["DRA", "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1"]

    classI_imgt_alleles = set()
    classII_imgt_alleles = set()
    for allele in missing_alleles:
        gene = allele.split("*")[0]
        if gene in class_I_genes:
            classI_imgt_alleles.add(allele)
        if gene in class_II_genes:
            classII_imgt_alleles.add(allele)

    mro_alleles = set()
    with open("../ontology/molecule.tsv") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            mro_protein = row["IEDB Label"]
            if "HLA" in mro_protein and "*" in mro_protein:
                if "/" in mro_protein:
                    mro_alleles.add(mro_protein.split("/")[0][3:])
                    mro_alleles.add(mro_protein.split("/")[1])
                else:
                    mro_alleles.add(mro_protein.split("HLA-")[1])

    missing_classI = classI_imgt_alleles.difference(mro_alleles)
    missing_classII = classII_imgt_alleles.difference(mro_alleles)
    
    new_molecules = set()
    for x in list(create_classI_prot(missing_classI)):
        new_molecules.add(x)
    for y in list(create_classII_prot(missing_classII)):
        new_molecules.add(y)
    
    return new_molecules


get_IMGT_data()
missing_alleles = update_chain_sequence()
missing_genes = update_chain(missing_alleles)
missing_loci = update_locus(missing_genes)
missing_molecules = update_molecules(missing_alleles)
update_index(missing_alleles, missing_genes, missing_loci, missing_molecules)

try:
    os.remove("hla_prot.fasta")
    os.remove("Allelelist.txt")
except:
    print("Unable to remove IMGT files")
