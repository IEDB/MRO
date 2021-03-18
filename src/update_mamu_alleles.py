import csv
import os
import sys
from subprocess import Popen, PIPE

def get_mamu_sequences():
    """Gets mamu sequences from IPD-MHC database"""
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
                if line.startswith(">IPD-MHC:NHP") and "Mamu" in line.split(" ")[1]:
                    accession = line.split(" ")[0][9:]
                    allele = line.split(" ")[1]
                    allele = (":").join(allele.split(":")[:2])
            else:
                seq += line.strip()
        seqs[accession] = seq
        allele_names[allele] = accession

    return seqs, allele_names

def get_current_loci():
    """Simple method to get current Mamu loci in MRO"""
    mro_loci = set()
    with open(sys.argv[4]) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            locus = row["Label"]
            if "Mamu" in locus:
                mro_loci.add(locus.split(" ")[0])

    return mro_loci

def update_chains(curr_loci, ipd_seqs, allele_map):
    """Updates the chain-sequence.tsv and chain.tsv with missing alleles from IPD

    Returns:
        A list of missing alleles ['A*02:01', 'B*02:01']
    """

    mro_alleles = set()
    with open(sys.argv[2]) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            if "Mamu" in row["Label"]:
                allele = row["Label"].split(" ")[0]
                allele = (":").join(allele.split(":")[:2])
                mro_alleles.add(allele)

    varied_loci = ["A", "DQA", "DPA", "DPB", "DRB", "DQB"]
    # Find differences between IMGT and MRO alleles
    imgt_alleles = set(allele_map.keys())
    missing_alleles = imgt_alleles.difference(mro_alleles)
    new_alleles = {x for x in missing_alleles if x.split("*")[0] in curr_loci or any(loci in x.split("*")[0] for loci in varied_loci) and "DMA" not in x.split("*")[0]}
    new_alleles = {x for x in new_alleles if x[-1] != "N"}

    missing_chain_seq_rows = set()
    missing_chain_rows = set()
    for allele in new_alleles:
        gene = allele.split("*")[0]
        if "A" == gene.split("-")[1][0]:
            gene = "Mamu-A"
        if "DPA" in gene:
            gene = "Mamu-DPA"
        if "DRB" in gene:
            gene = "Mamu-DRB"
        if "DPB" in gene:
            gene = "Mamu-DPB"
        if "DQA" in gene:
            gene = "Mamu-DQA"
        if "DQB" in gene:
            gene = "Mamu-DQB"
        try:
            accession = allele_map[allele]
            source = "IPD"
            label = f"{allele} chain"
            resource_name = allele
            seq = ipd_seqs[allele_map[allele]]
            missing_chain_seq_rows.add((label, resource_name, source, accession, seq))
            
            synonyms = ""
            class_type = "subclass"
            parent = f"{gene} chain"
            missing_chain_rows.add((label, synonyms, class_type, parent, "", ""))
        except Exception:
            print(f"Allele {allele} has no sequence in IPD")

    with open(sys.argv[2], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in missing_chain_rows:
            writer.writerow(tup)

    with open(sys.argv[1], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in missing_chain_seq_rows:
            writer.writerow(tup)

    return new_alleles

def create_classI_prot(missing_chainI):
    """Creates entries for MHC class I molecules

    Returns:
        A list of new entries for molecule.tsv
    """
    new_classI_molecules = set()
    with open(sys.argv[3], "a+") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for allele in missing_chainI:
            label = f"{allele} protein complex"
            iedb_name = f"{allele}"
            synonym = ""
            restrict_lvl = "complete molecule"
            class_type = "equivalent"
            parent = "MHC class I protein complex"
            taxon = "rhesus macaque"
            alpha_chain = f"{allele} chain"
            beta_chain = "Beta-2-microglobulin"
            writer.writerow(
                [
                    label,
                    iedb_name,
                    synonym,
                    restrict_lvl,
                    class_type,
                    parent,
                    taxon,
                    alpha_chain,
                    beta_chain,
                ]
            )
            new_classI_molecules.add(f"{allele} protein complex")
    return new_classI_molecules


def create_classII_pairing(allele):
    """Simple helper class to create a pairing of alpha and beta chain based on gene

    Returns:
        A string that has a pairing of alpha and beta chain (i.e. DRBA/DRB1*01:02)
    """
    pairing = {
        "Mamu-DPA1": "Mamu-DPB",
        "Mamu-DRA": "Mamu-DRB",
        "Mamu-DQA1": "Mamu-DQB"
    }
    pairing_rev = {
        "Mamu-DQB1": "Mamu-DQA",
        "Mamu-DPB1": "Mamu-DPA",
        "Mamu-DPB2": "Mamu-DPA",
        "Mamu-DRB": "Mamu-DRA",
        "Mamu-DRB1": "Mamu-DRA",
        "Mamu-DRB3": "Mamu-DRA",
        "Mamu-DRB4": "Mamu-DRA",
        "Mamu-DRB5": "Mamu-DRA",
        "Mamu-DRB6": "Mamu-DRA"
    }
    gene = allele.split("*")[0]
    if gene in pairing:
        pair = f"{allele}/{pairing[gene]}"
    else:
        pair = f"{pairing_rev[gene]}/{allele}"

    return pair


def create_classII_prot(missing_chainII):
    """Creates entries for MHC class II molecules

    Returns:
        A list of new entries for molecule.tsv
    """
    new_classII_molecules = set()
    with open(sys.argv[3], "a+") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for allele in missing_chainII:
            pair = create_classII_pairing(allele)
            alpha_gene = pair.split("/")[0]
            beta_gene = pair.split("/")[1]

            # partial molecule entry
            label = f"{allele} protein complex"
            iedb_name = f"{allele}"
            synonym = ""
            restrict_lvl = "partial molecule"
            class_type = "equivalent"
            parent = "MHC class II protein complex"
            taxon = "rhesus macaque"
            alpha_chain = f"{alpha_gene} chain"
            beta_chain = f"{beta_gene} chain"

            writer.writerow(
                [
                    label,
                    iedb_name,
                    synonym,
                    restrict_lvl,
                    class_type,
                    parent,
                    taxon,
                    alpha_chain,
                    beta_chain,
                ]
            )
            new_classII_molecules.add(f"{allele} protein complex")

    return new_classII_molecules

def update_molecules(missing_alleles):
    """Builds and updates Mamu molecules"""
    class_I_genes = ["A", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "AG", "B", "E"]
    class_II_genes = ["DPA1", "DPB1", "DRA", "DRB", "DRB1", "DRB3", "DRB4", "DRB5", "DRB6", "DQA1", "DQB1", "DPB2"]

    mro_proteins = set()
    with open(sys.argv[3]) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            mro_protein = row["IEDB Label"]
            if "Mamu" in mro_protein:
                if "/" in mro_protein:
                    continue
                else:
                    mro_proteins.add(mro_protein)
    
    classI = set()
    classII = set()
    for allele in missing_alleles:
        gene = allele.split("*")[0][5:]
        if gene in class_I_genes:
            classI.add(allele)
        if gene in class_II_genes:
            classII.add(allele)
    
    new_molecules = set()
    for x in list(create_classI_prot(classI)):
        if x not in mro_proteins:
            new_molecules.add(x)
    for y in list(create_classII_prot(classII)):
        if y not in mro_proteins:
            new_molecules.add(y)

    return new_molecules

def update_index(missing_alleles, missing_molecules):
    """Adds new entries from IPD to index.tsv to allow ROBOT to build owl file"""
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

    with open(sys.argv[5], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in new_tups:
            writer.writerow(tup)

def update_IEDB_tab(missing_molecules):
    """Increments IEDB IDs and adds new molecules"""
    # Get current IEDB sheet
    with open(sys.argv[7]) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
        final_row = rows[-1]
        curr_iedb_id = int(final_row["IEDB ID"])

    # Write out results
    with open(sys.argv[7], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for molecule in missing_molecules:
            locus = molecule.split("-")[1].split("*")[0]
            curr_iedb_id += 1
            # Specific to mamu alleles, either DR, DP or whatever locus is
            if "DR" in locus:
                tup = (molecule, curr_iedb_id, "DR", "", "")
            if "DP" in locus:
                tup = (molecule, curr_iedb_id, "DP", "", "")
            if "DQ" in locus:
                tup = (molecule, curr_iedb_id, "DQ", "", "")
            else:
                tup = (molecule, curr_iedb_id, locus, "", "")
            writer.writerow(tup)

ipd_seqs, allele_map = get_mamu_sequences()
curr_loci = get_current_loci()
missing_alleles = update_chains(curr_loci, ipd_seqs, allele_map)
missing_molecules = update_molecules(missing_alleles)
update_index(missing_alleles, missing_molecules)
update_IEDB_tab(missing_molecules)