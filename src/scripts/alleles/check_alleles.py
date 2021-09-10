import json


alleles = []
with open("~/MRO/ontology/genetic-locus.tsv") as f:
    for line in f:
        if "-" in line:
            alleles.append(line.split("\t")[0].split(" ")[0])

loci = []
with open("MRO_loci_to_add.txt") as f:
    for line in f:
        loci.append(line.rstrip())

missing_loci = []
for locus in loci:
    if locus not in alleles:
        missing_loci.append(locus)

with open("MRO_loci_to_add_unique.txt", "w") as f:
    for locus in missing_loci:
        f.write(locus + "\n")
