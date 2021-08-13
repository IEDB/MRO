import csv
import os
import sys
from subprocess import Popen, PIPE
import json


def get_species_metadata():
    with open("src/scripts/alleles/locus_data.json", "r") as fh:
        data = json.load(fh)

    return data


def get_sequences():
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
                    allele = (":").join(allele.split(":")[:2])
            else:
                seq += line.strip()
        seqs[accession] = seq
        allele_names[allele] = accession
    return seqs, allele_names


def get_taxon_metadata(metadata):
    # Map genes to taxon
    taxon_metadata = {}
    for species in metadata["species"]:
        mhc_code = species["mhc_code"]
        taxon = species["taxon"]
        for gene in species["class_I"]:
            taxon_metadata[mhc_code + "-" + gene] = taxon
        for gene in species["pairing"].keys():
            alpha = gene
            beta = species["pairing"][gene]
            taxon_metadata[mhc_code + "-" + alpha] = taxon
            taxon_metadata[mhc_code + "-" + beta] = taxon

    return taxon_metadata


def get_pairings(metadata):
    pairings = {}
    pairings_rev = {}
    for species in metadata["species"]:
        mhc_code = species["mhc_code"]
        for key in species["pairing"]:
            pairings[mhc_code + "-" + key] = mhc_code + "-" + species["pairing"][key]
        for key in species["pairing_rev"]:
            pairings_rev[mhc_code + "-" + key] = (
                mhc_code + "-" + species["pairing_rev"][key]
            )

    return pairings, pairings_rev


def get_current_loci():
    """Simple method to get current loci in MRO"""
    mro_loci = set()
    with open(sys.argv[4]) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            locus = row["Label"]
            mro_loci.add(locus.split(" ")[0])

    return mro_loci


def get_loci_classes():
    classI = set()
    classII = set()
    with open(sys.argv[4]) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            locus = row["Label"]
            type = row["Parent"]
            if "MHC class I locus" in type:
                classI.add(locus.split(" ")[0])
            if "MHC class II locus" in type:
                classII.add(locus.split(" ")[0])

    return classI, classII


def update_loci(current_loci, metadata):
    """Modifies the genetic_locus.tsv file with any new loci that are missing
    Iterates through JSON metadata to find missing loci, build rows for tsv using this data
    """
    new_loci_rows = set()
    current_loci = set(current_loci)
    missing_loci = set()

    for species in metadata["species"]:
        mhc_code = species["mhc_code"]
        taxon = species["taxon"]

        # Class I loci
        for classI in species["class_I"]:
            allele = mhc_code + "-" + classI
            if allele not in current_loci:
                missing_loci.add(allele)
                mhc_class = "MHC class I"
                label = f"{allele} locus"
                synonyms = ""
                class_type = "subclass"
                parent = f"{taxon} {mhc_class} locus"
                new_loci_rows.add((label, synonyms, class_type, parent, ""))

        # Class II loci
        for classII in species["pairing"].keys():
            allele1 = mhc_code + "-" + classII
            allele2 = mhc_code + "-" + species["pairing"][classII]
            if allele1 not in current_loci:
                missing_loci.add(allele1)
                mhc_class = "MHC class II"
                label = f"{allele1} locus"
                synonyms = ""
                class_type = "subclass"
                parent = f"{taxon} {mhc_class} locus"
                new_loci_rows.add((label, synonyms, class_type, parent, ""))
            if allele2 not in current_loci:
                missing_loci.add(allele2)
                mhc_class = "MHC class II"
                label = f"{allele2} locus"
                synonyms = ""
                class_type = "subclass"
                parent = f"{taxon} {mhc_class} locus"
                new_loci_rows.add((label, synonyms, class_type, parent, ""))

    with open(sys.argv[4], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in new_loci_rows:
            writer.writerow(tup)

    return missing_loci


def update_chains(curr_loci, ipd_seqs, allele_map, metadata, missing_loci):
    """Updates the chain-sequence.tsv and chain.tsv with missing alleles from IPD

    Returns:
        A list of missing alleles ['A*02:01', 'B*02:01']
    """

    mro_alleles = set()
    mro_gen_chains = set()
    with open(sys.argv[2]) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            allele = row["Label"].split(" ")[0]
            allele = (":").join(allele.split(":")[:2])
            mro_alleles.add(allele)
            if "*" not in row["Label"]:
                mro_gen_chains.add(row["Label"])

    # Load MHCFlurry alleles
    with open("src/scripts/alleles/mhcflurry_alleles.json", "r") as fh:
        data = json.load(fh)
    imgt_alleles = set(data)
    missing_alleles = imgt_alleles.difference(mro_alleles)
    new_alleles = {x for x in missing_alleles if x.split("*")[0] in curr_loci}
    new_alleles = {x for x in new_alleles if x[-1] != "N"}
    new_alleles = {x for x in new_alleles if "HLA" not in x}

    new_gen_chains = set()
    missing_chain_seq_rows = set()
    missing_chain_rows = set()
    missing_gen_chain_rows = set()
    failed_alleles = set()
    for allele in new_alleles:
        gene = allele.split("*")[0]
        if gene + " chain" not in mro_gen_chains:
            missing_gen_chain_rows.add(
                (gene + " chain", "", "equivalent", "protein", gene + " locus")
            )
            new_gen_chains.add(gene)
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
            failed_alleles.add(allele)
            print(f"Allele {allele} has no sequence in IPD")

    with open(sys.argv[2], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in missing_chain_rows:
            writer.writerow(tup)
        for tup in missing_gen_chain_rows:
            writer.writerow(tup)

    with open(sys.argv[1], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in missing_chain_seq_rows:
            writer.writerow(tup)

    new_alleles = new_alleles.difference(failed_alleles)

    return new_alleles, new_gen_chains


def create_classI_prot(missing_chainI):
    """Creates entries for MHC class I molecules

    Returns:
        A list of new entries for molecule.tsv
    """
    new_classI_molecules = set()
    with open(sys.argv[3], "a+") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for tup in missing_chainI:
            allele = tup[0]
            label = f"{allele} protein complex"
            iedb_name = f"{allele}"
            synonym = ""
            restrict_lvl = "complete molecule"
            class_type = "equivalent"
            parent = "MHC class I protein complex"
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
                    tup[1],
                    alpha_chain,
                    beta_chain,
                    "",
                    "",
                ]
            )
            new_classI_molecules.add(f"{allele} protein complex")
    return new_classI_molecules


def create_classII_pairing(allele, pairings, pairings_rev):
    """Simple helper class to create a pairing of alpha and beta chain based on gene

    Returns:
        A string that has a pairing of alpha and beta chain (i.e. DRBA/DRB1*01:02)
    """
    gene = allele.split("*")[0]
    if gene in pairings.keys():
        pair = f"{allele}/{pairings[gene]}"
    else:
        pair = f"{pairings_rev[gene]}/{allele}"

    return pair


def create_classII_prot(missing_chainII, pairings, pairings_rev):
    """Creates entries for MHC class II molecules

    Returns:
        A list of new entries for molecule.tsv
    """
    new_classII_molecules = set()
    with open(sys.argv[3], "a+") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for tup in missing_chainII:
            allele = tup[0]
            pair = create_classII_pairing(allele, pairings, pairings_rev)
            alpha_gene = pair.split("/")[0]
            beta_gene = pair.split("/")[1]

            # partial molecule entry
            label = f"{allele} protein complex"
            iedb_name = f"{allele}"
            synonym = ""
            restrict_lvl = "partial molecule"
            class_type = "equivalent"
            parent = "MHC class II protein complex"
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
                    tup[1],
                    alpha_chain,
                    beta_chain,
                    "",
                    "",
                ]
            )
            new_classII_molecules.add(f"{allele} protein complex")

    return new_classII_molecules


def create_classI_generic(mro_proteins, gene, taxon):
    if f"{gene} chain" not in mro_proteins:
        with open(sys.argv[3], "a+") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")

            label = f"{gene} protein complex"
            iedb_name = f"{gene}"
            synonym = ""
            restrict_lvl = "locus"
            class_type = "equivalent"
            parent = "MHC class I protein complex"
            alpha_chain = f"{gene} chain"
            beta_chain = f"Beta-2-microglobulin"

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
                    "",
                    "",
                ]
            )

    new_tups = []
    mro_labels = set()
    with open(sys.argv[5], "r") as fh:
        for _ in range(2):
            next(fh)
        for line in fh:
            curr_mro_id = line.rstrip().split("\t")[0]
            curr_mro_label = line.rstrip().split("\t")[1]
            mro_labels.add(curr_mro_label)
            curr_mro_num = int(curr_mro_id.split(":")[1])

        molecule_name = f"{gene} protein complex"
        if molecule_name not in mro_labels:
            mro_id = f"MRO:{str(curr_mro_num + 1).zfill(7)}"
            curr_mro_num += 1
            new_tups.append((mro_id, molecule_name, "owl:Class"))

    with open(sys.argv[5], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in new_tups:
            writer.writerow(tup)


def create_classII_generic(mro_proteins, gene, taxon):
    if f"{gene} chain" not in mro_proteins:
        with open(sys.argv[3], "a+") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")

            pair = create_classII_pairing(gene, pairings, pairings_rev)
            alpha_gene = pair.split("/")[0]
            beta_gene = pair.split("/")[1]

            # partial molecule entry
            label = f"{gene} protein complex"
            iedb_name = f"{gene}"
            synonym = ""
            restrict_lvl = "partial molecule"
            class_type = "equivalent"
            parent = "MHC class II protein complex"
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
                    "",
                    "",
                ]
            )

    new_tups = []
    mro_labels = set()
    with open(sys.argv[5], "r") as fh:
        for _ in range(2):
            next(fh)
        for line in fh:
            curr_mro_id = line.rstrip().split("\t")[0]
            curr_mro_label = line.rstrip().split("\t")[1]
            mro_labels.add(curr_mro_label)
            curr_mro_num = int(curr_mro_id.split(":")[1])

        molecule_name = f"{gene} protein complex"
        if molecule_name not in mro_labels:
            mro_id = f"MRO:{str(curr_mro_num + 1).zfill(7)}"
            curr_mro_num += 1
            new_tups.append((mro_id, molecule_name, "owl:Class"))

    with open(sys.argv[5], "a+") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for tup in new_tups:
            writer.writerow(tup)


def update_molecules(
    missing_alleles, classI, classII, taxon_metadata, pairings, pairings_rev
):
    """Builds and updates MHC molecules"""
    mro_proteins = set()
    with open(sys.argv[3]) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            mro_protein = row["IEDB Label"]
            if "/" in mro_protein:
                continue
            else:
                mro_proteins.add(mro_protein)

    class_I = set()
    class_II = set()
    for allele in missing_alleles:
        gene = allele.split("*")[0]
        if gene in classI:
            create_classI_generic(mro_proteins, gene, taxon_metadata[gene])
            class_I.add((allele, taxon_metadata[gene]))
        if gene in classII:
            create_classII_generic(mro_proteins, gene, taxon_metadata[gene])
            class_II.add((allele, taxon_metadata[gene]))

    new_molecules = set()
    for x in list(create_classI_prot(class_I)):
        if x not in mro_proteins:
            new_molecules.add(x)
    for y in list(create_classII_prot(class_II, pairings, pairings_rev)):
        if y not in mro_proteins:
            new_molecules.add(y)

    return new_molecules


def update_index(missing_alleles, missing_molecules, missing_loci, new_gen_chains):
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
            if "DQ" in locus:
                tup = (molecule, curr_iedb_id, "DQ", "", "")
            elif "DRB" in locus:
                tup = (molecule, curr_iedb_id, "DRB", "", "")
            else:
                tup = (molecule, curr_iedb_id, locus, "", "")
            writer.writerow(tup)


metadata = get_species_metadata()
pairings, pairings_rev = get_pairings(metadata)
taxon_metadata = get_taxon_metadata(metadata)
ipd_seqs, allele_map = get_sequences()
curr_loci = get_current_loci()
missing_loci = update_loci(curr_loci, metadata)
# get current loci again
curr_loci = get_current_loci()
classI, classII = get_loci_classes()
missing_alleles, new_gen_chains = update_chains(
    curr_loci, ipd_seqs, allele_map, metadata, missing_loci
)
missing_molecules = update_molecules(
    missing_alleles, classI, classII, taxon_metadata, pairings, pairings_rev
)
update_index(missing_alleles, missing_molecules, missing_loci, new_gen_chains)
update_IEDB_tab(missing_molecules)
