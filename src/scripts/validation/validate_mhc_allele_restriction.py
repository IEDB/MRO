import csv
import sys

from argparse import ArgumentParser
from cerberus import Validator

ALLOWED = "allowed"
MAX_LEN = "maxlength"
RE = "regex"
TYPE = "type"
NULL = "nullable"


def main():
    p = ArgumentParser()
    p.add_argument("mhc_allele_restriction")
    p.add_argument("output")
    args = p.parse_args()

    schema = {
        "mhc_allele_restriction_id": {TYPE: "number", NULL: False},
        "displayed_restriction": {TYPE: "string", MAX_LEN: 85, NULL: False},
        "synonyms": {TYPE: "string", MAX_LEN: 200, NULL: True},
        "includes": {TYPE: "string", MAX_LEN: 85, NULL: True},
        "restriction_level": {
            TYPE: "string",
            MAX_LEN: 35,
            ALLOWED: [
                "class",
                "complete molecule",
                "haplotype",
                "locus",
                "partial molecule",
                "serotype",
            ],
            NULL: False,
        },
        "organism": {TYPE: "string", MAX_LEN: 150, NULL: False},
        "organism_ncbi_tax_id": {TYPE: "number", NULL: False},
        "class": {
            TYPE: "string",
            MAX_LEN: 35,
            ALLOWED: ["I", "II", "non classical"],
            NULL: False,
        },
        "haplotype": {TYPE: "string", MAX_LEN: 10, NULL: True},
        "locus": {MAX_LEN: 10, NULL: True},
        "serotype": {TYPE: "string", MAX_LEN: 10, NULL: True},
        "molecule": {TYPE: "string", MAX_LEN: 50, NULL: True},
        "chain_i_name": {TYPE: "string", MAX_LEN: 35, NULL: True},
        "chain_ii_name": {TYPE: "string", MAX_LEN: 35, NULL: True},
        "chain_i_locus": {MAX_LEN: 10, NULL: True},
        "chain_i_mutation": {TYPE: "string", MAX_LEN: 35, NULL: True},
        "chain_ii_locus": {TYPE: "string", MAX_LEN: 10, NULL: True},
        "chain_ii_mutation": {TYPE: "string", MAX_LEN: 35, NULL: True},
        "chain_i_source_id": {TYPE: "number", NULL: True},
        "chain_ii_source_id": {TYPE: "number", NULL: True},
        "chain_i_accession": {TYPE: "string", RE: r"HLA[0-9]+", NULL: True},
        "chain_ii_accession": {TYPE: "string", RE: r"HLA[0-9]+", NULL: True},
        "chain_i_mro_id": {TYPE: "string", RE: r"MRO:[0-9]{7}", NULL: True},
        "chain_ii_mro_id": {TYPE: "string", RE: r"MRO:[0-9]{7}", NULL: True},
        "iri": {
            TYPE: "string",
            MAX_LEN: 100,
            RE: r"^http://purl\.obolibrary\.org/obo/MRO_[0-9]{7}$",
            NULL: False,
        },
    }
    validator = Validator(schema)
    errs = []
    line = 2
    with open(args.mhc_allele_restriction, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            fixed_row = {}
            for k, v in row.items():
                if v == "":
                    fixed = None
                else:
                    try:
                        fixed = int(v)
                    except Exception as e:
                        fixed = v
                fixed_row[k] = fixed
            if not validator.validate(fixed_row):
                for col_name, msg in validator.errors.items():
                    fmt = {
                        "Line": line,
                        "Column": col_name,
                        "Value": row[col_name],
                        "Message": ", ".join(msg),
                    }
                    errs.append(fmt)
            line += 1

    with open(args.output, "w") as f:
        if len(errs) > 0:
            print(
                f"{len(errs)} error(s) found in {args.mhc_allele_restriction} - see {args.output}!"
            )
            writer = csv.DictWriter(
                f,
                delimiter="\t",
                lineterminator="\n",
                fieldnames=["Line", "Column", "Value", "Message"],
            )
            writer.writeheader()
            writer.writerows(errs)
            sys.exit(1)
        else:
            print(f"No errors found in {args.mhc_allele_restriction}")
            f.write("")


if __name__ == "__main__":
    main()
