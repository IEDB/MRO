import Bio
from Bio import SeqIO
import csv
import sys
import os
sys.path.append(os.getcwd())
from io import StringIO
from Bio import BiopythonWarning


chain_sequence = open("ontology/chain-sequence.tsv", "r")

reader = csv.DictReader(chain_sequence, delimiter = "\t")
next(reader)
data = {row["Accession"]: row["Label"] for row in reader}


acc = list(data.keys())
excluded_sequence = []

G_domains = []
for entry in SeqIO.parse("build/hla.dat", "imgt" ):
    has_exon_1 = False
    if entry.name not in acc or entry.description.split(",")[0].endswith("N"):
        continue
    exons = []
    cds = ""
    exons_to_get = []
    type = -1
    if 'Human MHC Class I sequence' in entry.description and 'partial' in entry.description:
        exons_to_get = ['2', '3']
        type = 1
    elif 'Human MHC Class II sequence' in entry.description and 'partial' in entry.description:
        exons_to_get = ['2']
        type = 2
    elif 'Human MHC Class I sequence' in entry.description:
        exons_to_get = ['2', '3']
        type = 3
    else:
        exons_to_get = ['2']
        type = 4
    G_domain_nuc = ""
    G_domain = ""
    try:
        for feature in entry.features:
            if 'number' in feature.qualifiers and feature.qualifiers['number'][0] == '1' and feature.type == 'exon':
                has_exon_1 = True
                exons.append(feature)
            if 'number' in feature.qualifiers and (feature.qualifiers['number'][0] in exons_to_get) and feature.type == 'exon':
                exons.append(feature)
            if feature.type == 'CDS':
                cds = feature
            if len(exons) == 3:
                break
        codon_start = int(cds.qualifiers['codon_start'][0])
        if type == 3 and len(exons) == 3:
            exon_1_len = exons[0].location.end - exons[0].location.start
            last_full_codon_start = list(range(codon_start - 1, exon_1_len -1, 3))
            last_full_codon_start = last_full_codon_start[-1]
            #except IndexError:
            #    import pdb; pdb.set_trace()
            splice_loc = last_full_codon_start + 2 + 1
            splice = exons[0].extract(entry)[splice_loc: ]
            exon_2 = exons[1].extract(entry)
            G_domain_nuc = splice + exon_2
            exon_3 = exons[2].extract(entry)
            G_domain_nuc = G_domain_nuc + exon_3
        elif type == 4 and len(exons) == 2:
            exon_1_len = exons[0].location.end - exons[0].location.start
            last_full_codon_start = list(range(codon_start - 1, exon_1_len -1, 3))[-1]
            splice_loc = last_full_codon_start + 2 + 1
            splice = exons[0].extract(entry)[splice_loc: ]
            exon_2 = exons[1].extract(entry)
            G_domain_nuc = splice + exon_2
        elif type == 2 and len(exons) == 2 and has_exon_1:
            exon_1_len = exons[0].location.end - exons[0].location.start
            if exon_1_len == 1:
                splice_loc = 0
            else:
                last_full_codon_start = list(range(codon_start - 1, exon_1_len -1, 3))[-1]
                splice_loc = last_full_codon_start + 2 + 1
            splice = exons[0].extract(entry)[splice_loc: ]
            exon_2 = exons[1].extract(entry)
            G_domain_nuc = splice + exon_2
        elif type == 2 and len(exons) == 1 and not has_exon_1:
            G_domain = cds.qualifiers['translation'][0]
        elif type == 1 and len(exons) == 3 and has_exon_1:
            exon_1_len = exons[0].location.end - exons[0].location.start
            if exon_1_len == 1:
                splice_loc = 0
            else:
                last_full_codon_start = list(range(codon_start - 1, exon_1_len -1, 3))[-1]
                splice_loc = last_full_codon_start + 2 + 1
            splice = exons[0].extract(entry)[splice_loc: ]
            exon_2 = exons[1].extract(entry)
            G_domain_nuc = splice + exon_2
            exon_3 = exons[2].extract(entry)
            G_domain_nuc = G_domain_nuc + exon_3
        elif type == 1 and len(exons) == 2 and not has_exon_1:
            G_domain = cds.qualifiers['translation'][0]
        else:
            #import pdb; pdb.set_trace()
            excluded_sequence.append(entry.name)
            continue
        if not G_domain:
            G_domain = G_domain_nuc.translate()
        if 'seq' in dir(G_domain) and str(G_domain.seq) not in cds.qualifiers['translation'][0]:
            import pdb; pdb.set_trace()
            print(entry.name, G_domain.seq, cds.qualifiers['translation'][0], "not matched" )
        if 'seq' in dir(G_domain):
            G_domain =  str(G_domain.seq)
        else:
            G_domain = G_domain
        G_domains.append({"Label": data[entry.name], "minimal G domain sequence" : G_domain})
    except BiopythonWarning:
        print("BiopythonWarning")
    except AttributeError:
        if str(entry.seq) == 'X':
            excluded_sequence.append(entry)

print(excluded_sequence)

import csv
with open("ontology/G-domain-sequence.tsv", "w") as fh:
    writer = csv.DictWriter(fh, fieldnames = G_domains[0].keys(), delimiter = "\t")
    writer.writeheader()
    fh.write("LABEL\tA minimal G domain sequence")
    writer.writerows(G_domains)
