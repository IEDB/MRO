import csv
import os
from Bio import SeqIO
# gen_records = {}
# records_list = [seq_record for seq_record in SeqIO.parse(os.path.join("/Users/amody/IMGTHLA/fasta",  "hla_gen.fasta"), "fasta")]
# for record in records_list:
#     stuff = record.description.split(" ")
#     record.id = stuff[0][4:]
#     record.name = stuff[1]
# gen_records.update(SeqIO.to_dict(records_list, key_function = lambda rec: rec.name))

nuc_records = {}
rec = SeqIO.parse("build/hla.dat", "imgt")
for b in rec:
    exons = []
    for feature in b.features:
        if 'number' in feature.qualifiers:
            number = feature.qualifiers['number']
            if feature.type == 'exon' and (number == ['2'] or number == ['3'] ):
                exons.append(str(feature.extract(b.seq)))
    exons = "|".join(exons)
    nuc_records[b.description.split(",")[0]] = exons
gen_seq = []
G_groups = []
p = open("build/hla_nom_g.txt", "r")
c = p.read().splitlines()
for i in c:
    if i.startswith("#"):
        continue
    yt = i.split(";")
    locus = yt[0]
    alleles = yt[1].split("/")
    if len(yt) > 2 and len(yt[2]) > 0:
        G_grp = "HLA-" + locus + yt[2]
        nuc_seqs = {"G group": G_grp, "seqs": set()}
        for allele in alleles:
            nuc_seqs["seqs"].add(str(nuc_records["HLA-" + locus + allele]))
            g_dict = {}
            g_dict["G group"] = G_grp
            g_dict["MHC gene allele"] = "HLA-" + locus + allele
            #if locus+allele not in gen_records.keys():
                #gen_seq.append(g_dict)
                #continue
            #g_dict["genomic sequence"] = str(gen_records[locus + allele].seq)
            #gen_seq.append(g_dict)
        if len(nuc_seqs["seqs"]) == 1:
            nuc_seqs["seqs"] = next(elem for elem in nuc_seqs["seqs"])
            nuc_seqs["subclass"] = "G group"
            G_groups.append(nuc_seqs)
        else:
            nuc_seqs["seqs"] = str(nuc_seqs["seqs"])
        nuc_seqs["subclass"] = "G group"
        G_groups.append(nuc_seqs)

file_obj = open("ontology/gene-alleles.tsv", "w")
writer = csv.DictWriter(file_obj, fieldnames = ["MHC gene allele", "G group", "genomic sequence"], delimiter = "\t")
writer.writeheader()
writer.writerows(gen_seq)
file_obj.close()


file_obj = open("ontology/G-group.tsv", "w")
writer = csv.DictWriter(file_obj, fieldnames = ["G group", "seqs", "subclass"], delimiter = "\t")
writer.writeheader()
writer.writerows(G_groups)
file_obj.close()
