import csv
import os
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError as TError
import json
import argparse
import sys
import re
import pandas as pd
from itertools import compress
# gen_records = {}
# records_list = [seq_record for seq_record in SeqIO.parse(os.path.join("/Users/amody/IMGTHLA/fasta",  "hla_gen.fasta"), "fasta")]
# for record in records_list:
#     stuff = record.description.split(" ")
#     record.id = stuff[0][4:]
#     record.name = stuff[1]
# gen_records.update(SeqIO.to_dict(records_list, key_function = lambda rec: rec.name))

EXCLUDED_GENES = {
    "BTN3A",
    "CD1",
    "DOA",
    "DOB",
    "V",
    "P",
    "S",
    "Y",
    "MICA",
    "MICB",
    "MR1",
    "HFE",
    "J",
    "L",
    "TAP1",
    "DMB",
    "DMA",
    "T",
    "V",
    "W",
    "K",
    "TAP2",
    "U",
    "H",
    "DRB8",
    "DRB9",
    "DRB7",
    "DRB6",
    "DPB2",
    "DPA2",
    "DQA2",
    "DRB2"
}
def get_chains():
    se = open("ontology/chain-sequence.tsv", "r")
    next(se)
    rows = csv.DictReader(se, delimiter="\t")
    chains = {row["LABEL"].split(" ")[0]: row["A MRO:sequence"] for row in rows}
    return chains
#nuc_records = {}

# m = []
# for i in rec:
#
#     types = [feature.type for feature in i.features]
#     if "CDS" not in types:
#         print(i.description)
#         m.append(str((i.name, i.description)))
# print("\n".join(m))

def get_G_groups():
    gen_seq = {}
    p = open("build/hla_nom_g.txt", "r")
    c = p.read().splitlines()
    for i in c:
        if i.startswith("#"):
            continue
        yt = i.split(";")
        locus = yt[0]
        alleles = yt[1].split("/")
        alleles = ["HLA-" + locus + allele for allele in alleles]
        if len(yt) > 2 and len(yt[2]) > 0:
            G_grp = "HLA-" + locus + yt[2]
            gen_seq.update({allele : G_grp for allele in alleles})
    return gen_seq

def update_allele_dict(allele_dict):
    allele_dict["MHC gene allele"] = allele_dict["MHC gene allele"] + " gene allele"
    allele_dict["Subclass"] = "MHC gene allele"
    return allele_dict

def process_hla_dat(gen_seq, gene_allele_fields,chains ):
    errors = []
    gen_alleles = []
    G_groups = {}
    rec = SeqIO.parse("build/hla.dat", "imgt")
    descriptions = set()
    for b in rec:
        try:
            gene_allele = {}
            desc = b.description.split(",")
            allele = desc[0]
            mhc_class = desc[1]
            mhc_class = ("II" if "II" in mhc_class else "I")
            gene_allele["MHC gene allele"] = allele
            elems = allele.split("*")
            locus = elems[0]
            if locus in EXCLUDED_GENES or locus.split("-")[1] in EXCLUDED_GENES or allele[len(allele) -1] == "N" or allele[len(allele) -1] == "Q":
                continue
            descriptions.add(desc[1])
            if mhc_class == "I":
                exons = [str(feature.extract(b).seq) for feature in b.features if feature.type == 'exon' and (feature.qualifiers['number'] == ['2'] or feature.qualifiers['number'] == ['3'])]
                exons = "|".join(exons)
            else:
                exons = [str(feature.extract(b).seq) for feature in b.features if feature.type == 'exon' and (feature.qualifiers['number'] == ['2'])]
                exons = exons[0]
            if allele in gen_seq and gen_seq[allele]:
                G_groups.setdefault(gen_seq[allele], set()).add(exons)
            cds = [feature for feature in b.features if feature.type=='CDS' and feature.location is not None and 'translation' in feature.qualifiers]
            if len(cds) == 1:
                cds = cds[0]
            else:
                continue
            gene_allele["Coding Region Sequence"] = str(cds.extract(b).seq)
            gene_allele["Source"] = "IMGT/HLA"
            gene_allele["Accession"] = b.name
            match = re.search(pattern = r"[0-9]+$", string=locus)
            if match:
                locus = locus[:match.span()[0]]

            gene_allele["Locus"] = locus + " locus"
            two_field = allele.split(":")
            two_field = two_field[0] + ":" + two_field[1]
            if two_field in chains:
                #if chains[two_field] == cds.extract(b).seq[int(cds.qualifiers['codon_start'][0])-1:].translate():
                if str(cds.translate(b).seq) in chains[two_field] or chains[two_field] in str(cds.translate(b).seq):
                    gene_allele["Chain"] = two_field + " chain"
                else:
                    error = {"reason":"MRO protein not equal to IMGT protein", "IMGT Accession": b.name, "MRO allele": two_field, "MRO protein": chains[two_field], "Translated Protein": str(cds.translate(b).seq)   }
                    errors.append(error)
            else:
                error = {"reason": "MRO doesn't have allele", "IMGT Accession": b.name, "IMGT allele" : allele, "MRO allele": two_field}
                errors.append(error)
        except AttributeError:
            error = {"reason" : "AttributeError", "IMGT Accession": b.name, "IMGT allele" : allele}
            errors.append(error)
            continue
        except IndexError:
            error = {"reason" : "IndexError", "IMGT Accession": b.name, "IMGT allele" : allele}
            errors.append(error)
            continue
        except TError:
            # mainly for alleles with partial sequences
            #print("TranslationError", b.name)
            extracted_protein = cds.qualifiers["translation"][0]
            extracted_protein = str(extracted_protein) 
            if chains[two_field] in extracted_protein or extracted_protein in chains[two_field]:
                gene_allele["Chain"] = two_field + " chain"
            else:
                error = {"reason" : "TranslationError", "IMGT Accession": b.name, "IMGT allele" : allele}
                errors.append(error)
        if allele in gen_seq:
            gene_allele["G group"] = f"'{gen_seq[allele]}'"
        else:
            gene_allele["G group"] = ''
        all_fields_present = True
        excluded_fields = ["G group", "MHC gene allele", "Subclass"]
        for field in gene_allele_fields:
            if field not in excluded_fields and field not in gene_allele.keys():
                all_fields_present = False
        if all_fields_present:
            gen_alleles.append(gene_allele)
        continue
    G_groups = [{"G group" : allele, "Exon 2 and/or 3": max(G_groups[allele], key=len), "Logic": "G group"} for allele in G_groups]
    gen_alleles = list(map(update_allele_dict, gen_alleles))
    with open("descriptions.txt", "w") as desc_file:
        for i in descriptions:
            desc_file.write(i + "\n")
    return gen_alleles, G_groups, errors

def write_error_report(errors):
    with open("build/report-g-grp.json", "w") as report:
        json.dump(errors, report)

def write_gene_alleles(gene_allele_fields, gen_alleles):
    with open("ontology/gene-alleles.tsv", "w") as file_obj:

        writer = csv.DictWriter(file_obj, fieldnames = gene_allele_fields, delimiter = "\t")
        writer.writeheader()
        #file_obj.write("LABEL\tEC 'has gene product' some %\tEC 'has part' some %\tA MRO:accession\tA MRO:source\tA MRO:sequence\n")
        #writer.writerows([gen_alleles[0]])
        file_obj.write("LABEL\tSC 'has gene product' some %\tSC 'has part' some %\tSC %\tSC 'gene product of' some %\tA MRO:accession\tA MRO:source\tA MRO:sequence\n")
        writer.writerows(gen_alleles)
        file_obj.close()

def write_G_groups(G_groups):
    with open("ontology/G-group.tsv", "w") as file_obj:
        writer = csv.DictWriter(file_obj, fieldnames = ["G group", "Exon 2 and/or 3", "Logic"], delimiter = "\t")
        writer.writeheader()
        file_obj.write("LABEL\tA MRO:sequence\tSC %\n")
        writer.writerows([G_groups[0]])
        writer.writerows(G_groups)
        file_obj.close()

def get_index_info():
    with open("index.tsv", "r") as index:
        entries = index.read().splitlines()
        index.close()
        entries = list(map(lambda x: x.rstrip(), entries))
        first = entries[0]
        second = entries[1]
        entries = entries[2:]
        entries = list(map(lambda entry: dict(zip(["ID", "Label", "Type", "Depreciated?"], entry.split("\t"))), entries))
        for entry in entries:
            print(entry)
            int(entry["ID"].split(":")[1])
        ids = [int(entry["ID"].split(":")[1]) for entry in entries]
        cur_mro_id = max(ids) + 1
        labels = [entry["Label"] for entry in entries]
    return entries, ids, cur_mro_id, labels
    
def check_pop_properties_are_in_index(pop_properties):
    properties = [[prop, "owl:AnnotationProperty", ""] for prop in pop_properties]
    return properties 
    
def get_new_alleles_and_G_groups(gen_alleles, G_groups):
    new_alleles = [[allele["MHC gene allele"], "owl:Class", ""] for allele in gen_alleles]
    new_alleles.append(["MHC gene allele", "owl:Class", ""])
    new_G_groups = [[allele["G group"], "owl:Class", ""]  for allele in G_groups]
    new_G_groups.append(["G group", "owl:Class", ""])
    return new_alleles, new_G_groups


def update_index(terms):
    entries, ids, cur_mro_id, labels = get_index_info()
    print("start", cur_mro_id)
    x = len(terms) - 1
    while x >= 0:
        if terms[x][0] not in labels:
            terms[x].insert(0, "MRO:" + str(cur_mro_id).zfill(7))
            cur_mro_id +=1
        else:
            terms.remove(terms[x])
        x = x - 1
    print("end", cur_mro_id)
    with open("index.tsv", "a+") as index:
        writer = csv.writer(index, delimiter="\t", lineterminator="\n")
        writer.writerows(terms)

def read_template_data():
    alleles_file = open("ontology/gene-alleles.tsv", "r")
    fields = next(alleles_file).replace("\n", "").split("\t")
    next(alleles_file)
    gene_alleles = csv.DictReader(alleles_file, fieldnames = fields, delimiter="\t")
    gene_alleles = [dict(row) for row in gene_alleles]
    return gene_alleles

def verify_accession_data(data, gene_alleles):
    primary, secondary = get_allele_id_label(data)
    m = pd.DataFrame(data.loc[:,  (primary, secondary)].dropna(), copy = True)
    m.columns = m.columns.droplevel(0)
    m = m.rename(columns = {secondary: "Accession"})
    correction = m.isin(gene_alleles)
    missed_alleles = m.loc[~correction.Accession, :]
    missed_alleles = pd.Series(missed_alleles.Accession, copy = True)
    missed_alleles = set(missed_alleles)
    with open("build/report-g-grp.json") as report:
        x = json.load(report)
        excluded_alleles = [allele["IMGT Accession"] for allele in x]
    missed_alleles.difference_update(excluded_alleles)
        
    return set(missed_alleles)
    
def lookup_imgt(missed_alleles):
    import dbfetch as dbf
    batches= []
    missed_alleles = list(missed_alleles)
    for x in range(0, len(missed_alleles), 200):
        if x+200 <= len(missed_alleles):
            batches.append(dbf.fetchBatch(db = "imgthla",idListStr= ",".join(missed_alleles[x: x+200])))
        else:
            batches.append(dbf.fetchBatch(db = "imgthla", idListStr = ",".join(missed_alleles[x:])))
    batches = "".join(batches)
    return batches

def check_missed_alleles(imgt_data):
    from io import StringIO
    handler = StringIO(imgt_data)
    missed_alleles = SeqIO.parse(handler, "imgt")
    missed_alleles = {b.name: b for b in missed_alleles}
    same_alleles = []
    diff_alleles = []
    past_imgt = SeqIO.parse("build/hla1.dat", "imgt")
    for allele in past_imgt:
        if allele.name in missed_alleles:
            cds = [feature for feature in allele.features if feature.type=='CDS' and feature.location is not None and 'translation' in feature.qualifiers]
            if len(cds) == 1:
                cds = cds[0]
                cds1 = [feature for feature in missed_alleles[allele.name].features if feature.type=='CDS' and feature.location is not None and 'translation' in feature.qualifiers]
                if len(cds1)==1 and str(cds.qualifiers["translation"][0]) == str(cds1[0].qualifiers['translation'][0]):
                    same_alleles.append(allele.name)
                else:
                    diff_alleles.append(allele.name)
            else:
                diff_alleles.append(allele.name)
                
        else:
            continue
    return same_alleles, diff_alleles
    
def get_allele_id_label(data):
    levels = list(map(list, zip(*data.columns)))
    primary = ""
    secondary = ""
    for i in levels[0]:
        if "Allele Summary" in i:
            primary = i
            break
    for i in levels[1]:
        if "AlleleID" in i:
            secondary = i
            break
        if "Allele ID" in i:
            secondary = i
            break
    return primary, secondary
def verify_G_groups(data, gene_alleles):
    primary, secondary= get_allele_id_label(data)
    selection = [(primary, secondary), (primary, "G group")]
    m = pd.DataFrame(data.dropna(subset = selection ).loc[:, selection], copy = True)
    m.columns = m.columns.droplevel(0)
    m = m.rename(columns = {secondary: "Accession"})
    m.set_index("Accession", inplace = True)
    m.columns.name = ''
    temp = pd.DataFrame(gene_alleles.set_index("Accession").loc[:, "G group"].str.replace("'", ""), columns = ["G group"])
    temp.columns.name = ''
    correction = m.isin(temp)
    if correction.loc[correction["G group"] == False].empty:
        return set()
    else: 
        missed_alleles = m.loc[~correction["G group"], :]
        missed_alleles = pd.Series(missed_alleles.index, copy = True)
        missed_alleles = set(missed_alleles)
        
        with open("build/report-g-grp.json") as report:
            x = json.load(report)
            excluded_alleles = [allele["IMGT Accession"] for allele in x]
        missed_alleles.difference_update(excluded_alleles)
        
        return missed_alleles
def get_frequency_label(data):
    primary = ""
    columns = list(data.columns.levels[0])
    for i in columns:
        if "Frequency" in i:
            primary = i 
            break
    return primary
def add_totals(data, pop_group_map):
    total = pd.DataFrame(data.loc[data.index.str.endswith("total"), get_frequency_label(data)], copy = True)
    
    foo = total.index.str.replace(" total", "")
    total.index.name = "LABEL"
    total.rename(columns = pop_group_map, index = dict(zip(total.index, foo)) , inplace = True)
    total = total.add_prefix("AT '").add_suffix("'^^xsd:float")
    chains = pd.DataFrame(total.loc[~total.index.str.endswith("G") | ~total.index.str.endswith("N")], copy = True)
    chains.index = "HLA-" + chains.index + " chain"
    G_group = pd.DataFrame(total.loc[total.index.str.endswith("G")], copy = True)
    G_group.index = "HLA-" + G_group.index
    entries, ids, cur_mro_id, labels = get_index_info()
    in_mro = chains.index.isin(labels)
    not_in_mro_chains = chains[~in_mro]
    chains = chains.loc[in_mro]
    in_mro = G_group.index.isin(labels)
    print(in_mro)
    not_in_mro_G_group= G_group[~in_mro]
    G_group = G_group.loc[in_mro]
    print(len(not_in_mro_chains))
    print(len(not_in_mro_G_group))
    return chains, G_group
def main():
    parser = argparse.ArgumentParser(description='Update MHC gene allele sequences and G groups or add frequency data')
    parser.add_argument("-u","--update", action='store_true', help = "Update G groups and coding region genomic sequences of HLA alleles from IMGT")
    parser.add_argument("-f", "--frequency", action = 'store_true', help = "Update the frequency of each HLA allele in IMGT in population groups with data from CIWD 3.0")
    args = parser.parse_args()
    if args.update:
        gen_seq = get_G_groups()
        chains = get_chains()
        gene_allele_fields = ["MHC gene allele", "Chain", "G group", "Subclass", "Locus","Accession","Source","Coding Region Sequence"]
        gen_alleles, G_groups, errors = process_hla_dat(gen_seq = gen_seq, gene_allele_fields = gene_allele_fields, chains = chains)
        write_error_report(errors)
        write_gene_alleles(gene_allele_fields = gene_allele_fields, gen_alleles = gen_alleles)
        write_G_groups(G_groups)
        new_alleles, new_G_groups = get_new_alleles_and_G_groups(gen_alleles = gen_alleles, G_groups = G_groups)
        update_index(terms = new_alleles)
        update_index(terms= new_G_groups)
    if args.frequency:
        try:
            import pandas as pd
            import glob
            datafiles = glob.glob("build/HLA-*-frequency.xlsx")
            gene_alleles = read_template_data()
            gene_alleles = pd.DataFrame(gene_alleles)
            gene_alleles.loc[:, "MHC gene allele"] = gene_alleles.loc[:, "MHC gene allele"].str.replace(" gene allele", "").str.replace("HLA-", "")
            gene_alleles.loc[:, "G group"] = gene_alleles.loc[:, "G group"].str.replace("[']HLA-[A-Z]*[0-9]*[*]", "", regex = True).str.replace("'", "")
            gene_alleles = gene_alleles.set_index("MHC gene allele")
            missed_alleles = set()
            chains = []
            G_groups = []
            pop_group_map = {"AFA": "African frequency", "API": "Asian Pacific Islander frequency", "EURO": "European frequency", "MENA": "Middle Eastern and North African frequency", "HIS": "Hispanic frequency", "NAM": "Native American frequency", "UNK": "Unknown group frequency", "Total": "Total frequency"  }
            prop = check_pop_properties_are_in_index(list(pop_group_map.values()) + ["CIWD population frequency"])
            update_index(prop)
            template_string = "\t".join(pop_group_map.values())
            with open("ontology/chain-frequencies.tsv", "w") as totals:
                totals.write("Chains\t")
                totals.write(template_string)
                totals.write("\n")
            with open("ontology/G-group-frequencies.tsv", "w") as totals:
                totals.write("G group\t")
                totals.write(template_string)
                totals.write("\n")
            for datafile in datafiles:
                print(datafile)
                data = pd.read_excel(io = datafile, header = [0, 1], index_col = 0)
                data = data.loc[data.index.dropna()]
                data.drop(index = data.loc[data.index.str.endswith("N") | data.index.str.endswith("Q")].index, inplace=True)
                missing_alleles =set()
                missing_alleles |= verify_accession_data(data, gene_alleles)
                missing_alleles |= verify_G_groups(data, gene_alleles)
                missed_alleles |= missing_alleles
                chain, G_group = add_totals(data, pop_group_map)
                chains.append(chain)
                print(G_group)
                G_groups.append(G_group)
            chains = pd.concat(chains)
            G_groups = pd.concat(G_groups)
            print(G_groups)
            chains.to_csv("ontology/chain-frequencies.tsv", sep="\t", mode='a+')
            G_groups.to_csv("ontology/G-group-frequencies.tsv", sep="\t", mode='a+')
            if missed_alleles:
                imgt_data = lookup_imgt(missed_alleles)
                if check_missed_alleles(imgt_data):
                    print("foo")
            print(missed_alleles, len(missed_alleles))
        except ModuleNotFoundError:
            print("Please install pandas")


if __name__ == "__main__":
    main()


# file_obj = open("G-group.tsv", "w")
# writer = csv.DictWriter(file_obj, fieldnames = ["G group", "seqs", "subclass"], delimiter = "\t")
# writer.writeheader()
# writer.writerows(G_groups)
# file_obj.close()
