import csv
import os
from subprocess import Popen, PIPE


def get_IMGT_data():
    imgt_prot_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_prot.fasta"
    imgt_allele_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt"
    process = Popen(['wget', imgt_prot_URL, imgt_allele_URL], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    print(stdout, stderr)

def update_chain_sequence():
    """Update the chain-sequence.tsv table with any missing alleles and their sequences
    Returns a list of missing alleles"""

    #Parse the current MRO alleles
    mro_alleles = set()
    with open("../ontology/chain-sequence.tsv") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            if "HLA" in row["Label"]:
                mro_alleles.add(row["Label"].split(" ")[0][4:])

    # Create a dictionary mapping the accession to sequence from IMGT
    seqs = {}
    with open("./hla_prot.fasta") as fasta:
        accession = None
        seq = ''
        for line in fasta:
            if line.startswith('>'):
                if accession:
                    seqs[accession] = seq
                accession = None
                seq = ''

                # Match the accession
                if line.startswith('>HLA:HLA'):
                    accession = line[5:13]
                else:
                    print("Bad accession:", line)
            else:
                seq += line.strip()
        seqs[accession] = seq

    imgt_dict = {}

    # Get the mapping between allele name and IMGT accession
    with open("./Allelelist.txt") as fh:
        #Skip IMGT header
        for i in range(6):
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

    # Create a set for the new rows to avoid redundancy
    missing_allele_rows = set()
    for allele in missing_alleles:
        try:
            accession = imgt_dict[allele]
            source = "IMGT/HLA"
            label = "HLA-" + allele + " chain"
            resource_name = allele
            seq = seqs[imgt_dict[allele]]
            missing_allele_rows.add((label, resource_name, source, accession, seq))
        except:
            continue
    
    # File will be unsorted and needs to be passed to sort.py
    with open("../ontology/chain-sequence.tsv", "a+") as outfile:
        for tup in missing_allele_rows:
            outfile.write(("\t").join(tup)  + "\n")

    return missing_alleles

def update_chain(missing_alleles):
    """Update the chain.tsv table to add new chains"""

    #Get genes from the IMGT alleles not in MRO
    missing_genes = set()
    for allele in missing_alleles:
        gene = allele.split("*")[0]
        missing_genes.add(gene)

    #Parse out the current MRO alleles and genes in chain.tsv
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
    #Create rows for missing genes
    new_gene_tups = set()
    for gene in new_genes:
        label = "HLA-" + gene + " chain"
        synonyms = ""
        class_type = "equivalent"
        parent = "protein"
        gene = "HLA-" + gene + " locus"
        new_gene_tups.add((label, synonyms, class_type, parent, gene))

    new_chain_tups = set()
    #Create rows for missing chains
    for allele in missing_alleles:
        label = "HLA-" + allele + " chain"
        synonyms = ""
        class_type = "subclass"
        #This removes the specific allele #
        parent = "HLA-" + allele.split("*")[0] + " chain"
        gene = ""
        new_chain_tups.add((label, synonyms, class_type, parent, gene))

    with open("../ontology/chain.tsv", "a+") as fh:
        for tup in new_gene_tups:
            fh.write(("\t").join(tup) + "\n")
        for tup in new_chain_tups:
            fh.write(("\t").join(tup) + "\n")

    return new_genes

def update_locus(missing_genes):
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
        # Fast and dirty way to check if classI, classII, etc
        mhc_class = ""
        for prefix in class_I:
            if prefix in locus:
                mhc_class = "MHC class I"
        for prefix in class_II:
            if prefix in locus:
                mhc_class = "MHC class II"
        if mhc_class == "":
            mhc_class = "non-classical MHC"

        label = "HLA-" + locus + " locus"
        synonyms = ""
        class_type = "subclass"
        parent = "human " + mhc_class + " locus"
        new_loci_rows.add((label, synonyms, class_type, parent))
    
    with open("../ontology/genetic-locus.tsv", "a+") as fh:
        for row in new_loci_rows:
            fh.write(("\t").join(row)  + "\n")

    return new_loci

def update_index(missing_alleles, missing_genes, missing_loci):
    with open("../index.tsv") as fh:
        last_line = fh.readlines()[-1]
        curr_mro_label = last_line.rstrip().split("\t")[0]
        curr_mro_num = int(curr_mro_label.split(":")[1])
    
    new_tups = []
    for allele in missing_alleles:
        mro_id = "MRO:" + str(curr_mro_num + 1).zfill(7)
        curr_mro_num += 1
        chain_name = "HLA-" + allele + " chain"
        new_tups.append((mro_id, chain_name, "owl:Class"))

    for gene in missing_genes:
        mro_id = "MRO:" + str(curr_mro_num + 1).zfill(7)
        curr_mro_num += 1
        gene_name = "HLA-" + gene + " chain"
        new_tups.append((mro_id, gene_name, "owl:Class"))

    for locus in missing_loci:
        mro_id = "MRO:" + str(curr_mro_num + 1).zfill(7)
        curr_mro_num += 1
        locus_name = "HLA-" + locus + " locus"
        new_tups.append((mro_id, locus_name, "owl:Class"))

    with open("../index.tsv", "a+") as fh:
        for tup in new_tups:
            fh.write(("\t").join(tup) + "\n")
    

get_IMGT_data()
missing_alleles = update_chain_sequence()
missing_genes = update_chain(missing_alleles)
missing_loci = update_locus(missing_genes)
update_index(missing_alleles, missing_genes, missing_loci)

try:
    os.remove("hla_prot.fasta")
    os.remove("Allelelist.txt")
except:
    print("Unable to remove IMGT files")