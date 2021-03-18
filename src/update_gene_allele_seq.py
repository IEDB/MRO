import csv
import os
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError as TError
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
gene_allele_fields = ["MHC gene allele", "chain", "G group", "sequence"]
se = open("ontology/chain-sequence.tsv", "r")
next(se)
rows = csv.DictReader(se, delimiter="\t")
chains = {row["LABEL"].split(" ")[0]: row["A MRO:sequence"] for row in rows}
#nuc_records = {}
gen_alleles = []
gen_seq = {}
G_groups = {}
rec = SeqIO.parse("build/hla.dat", "imgt")
# m = []
# for i in rec:
#
#     types = [feature.type for feature in i.features]
#     if "CDS" not in types:
#         print(i.description)
#         m.append(str((i.name, i.description)))
# print("\n".join(m))


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
    else:
        gen_seq.update({allele : "" for allele in alleles})
for b in rec:
    try:
        gene_allele = {}
        allele, mhc_class = b.description.split(",")
        mhc_class = ("II" if "II" in mhc_class else "I")
        gene_allele["MHC gene allele"] = allele
        if allele.split("*")[0].split("-")[1] in EXCLUDED_GENES or allele.endswith("N"):
            continue
        if mhc_class == "I":
            exons = [str(feature.extract(b).seq) for feature in b.features if feature.type == 'exon' and (feature.qualifiers['number'] == ['2'] or feature.qualifiers['number'] == ['3'])]
            exons = "|".join(exons)
        else:
            exons = [str(feature.extract(b).seq) for feature in b.features if feature.type == 'exon' and (feature.qualifiers['number'] == ['2'])]
            
        if allele in gen_seq and gen_seq[allele]:
            G_groups.setdefault(gen_seq[allele], set()).add(exons)
        cds = [feature for feature in b.features if feature.type=='CDS' and feature.location is not None and 'translation' in feature.qualifiers]
        if len(cds) == 1:
            cds = cds[0]
        else:
            continue
        gene_allele["sequence"] = str(cds.extract(b).seq)
        two_field = allele.split(":")
        two_field = two_field[0] + ":" + two_field[1]
        if two_field in chains:
            #if chains[two_field] == cds.extract(b).seq[int(cds.qualifiers['codon_start'][0])-1:].translate():
            if str(cds.translate(b).seq) in chains[two_field] or chains[two_field] in str(cds.translate(b).seq):
                gene_allele["chain"] = two_field
            else:
                print("not equal protein", b.name, chains[two_field], two_field, str(cds.translate(b).seq) )
        else:
            print("doesn't have allele", b.name, allele, two_field)
    except AttributeError:
        print("AttributeError", b.name, allele)
        continue
    except IndexError:
        print("IndexError", b.name, allele)
        continue
    except TError:
        #print("TranslationError", b.name)
        if chains[two_field] == cds.extract(b).seq[int(cds.qualifiers['codon_start'][0])-1:].translate():
            gene_allele["chain"] = two_field
    if allele in gen_seq:
        gene_allele["G group"] = gen_seq[allele]
    else:
        gene_allele["G group"] = '""'
    all_fields_present = True
    for field in gene_allele_fields:
        if field not in gene_allele.keys():
            all_fields_present = False
    if all_fields_present:
        gen_alleles.append(gene_allele)
    continue

G_groups = [{"G group" : allele, "Exon 2 and 3": max(G_groups[allele], key=len)} for allele in G_groups]
with open("ontology/gene-alleles.tsv", "w") as file_obj:

    writer = csv.DictWriter(file_obj, fieldnames = gene_allele_fields, delimiter = "\t")
    writer.writeheader()
    file_obj.write("LABEL\tC 'encodes' some %\tC 'has part' some %\tA MRO:sequence\n")
    #writer.writerows([gen_alleles[0]])
    writer.writerows(gen_alleles)
    file_obj.close()
def update_g_dict(g_dict):
    g_dict.update({"Logic": "G group" })
    return g_dict
with open("ontology/G-group.tsv", "w") as file_obj:
    writer = csv.DictWriter(file_obj, fieldnames = ["G group", "Exon 2 and 3", "Logic"], delimiter = "\t")
    G_groups = list(map(lambda g_dict: update_g_dict(g_dict), G_groups))
    writer.writeheader()
    file_obj.write("LABEL\tA MRO:sequence\tEC %\n")
    print(G_groups[1])
    writer.writerows([G_groups[0]])
    writer.writerows(G_groups)
    file_obj.close()

# file_obj = open("G-group.tsv", "w")
# writer = csv.DictWriter(file_obj, fieldnames = ["G group", "seqs", "subclass"], delimiter = "\t")
# writer.writeheader()
# writer.writerows(G_groups)
# file_obj.close()
