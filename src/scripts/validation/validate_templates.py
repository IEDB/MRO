import csv
import re
import sys

from argparse import ArgumentParser


err_id = 0


def a1_to_idx(a1):
    """Convert an A1 cell label to row & column indexes. Adapted from gspread.utils."""
    m = re.compile(r"([A-Za-z]+)([1-9]\d*)").match(a1)
    if m:
        column_label = m.group(1).upper()
        row = int(m.group(2))

        col = 0
        for i, c in enumerate(reversed(column_label)):
            col += (ord(c) - 64) * (26 ** i)
    else:
        raise Exception(f"ERROR: Unable to parse cell label '{a1}'")

    return row, col


def check_labels(
    table_name,
    reader,
    label_source,
    valid_labels,
    regex=None,
    missing_level="error",
    allow_missing=False
):
    """Check that the labels used in a table are all present in a set of valid labels. If provided,
    also match the labels to a regex pattern. Return all labels from the table and a set of errors,
    if any."""
    global err_id
    errors = []
    row_idx = 3
    headers = reader.fieldnames
    labels = []

    missing_message = f"use a label defined in {label_source}"

    for row in reader:
        label = row["Label"]
        labels.append(label)
        if not allow_missing and label not in valid_labels:
            err_id += 1
            errors.append(
                {
                    "table": table_name,
                    "cell": idx_to_a1(row_idx, headers.index("Label") + 1),
                    "level": missing_level,
                    "rule ID": "unknown_label",
                    "rule": "unknown label",
                    "message": missing_message,
                }
            )
        if regex and not re.match(regex, label):
            err_id += 1
            errors.append(
                {
                    "table": table_name,
                    "cell": idx_to_a1(row_idx, headers.index("Label") + 1),
                    "level": "error",
                    "rule ID": "invalid_label",
                    "rule": "invalid label",
                    "message": f"change label to match pattern '{regex}'",
                }
            )
        row_idx += 1
    return labels, errors


def check_fields(
    table_name,
    reader,
    valid_labels,
    field_name="Parent",
    top_terms=None,
    source=None,
    required=True,
):
    """Validate that the contents of a given field (default=Parent) are present in a set of valid
    labels. If required, validate that this field value is filled in for all rows. Return a set of
    errors, if any."""
    global err_id
    if not top_terms:
        top_terms = []
    errors = []
    row_idx = 3
    headers = reader.fieldnames
    if not source:
        source = table_name
    for row in reader:
        value = row[field_name]
        if value and value.strip() == "":
            value = None

        if value:
            value = value.strip()

        if required and not value:
            err_id += 1
            errors.append(
                {
                    "table": table_name,
                    "cell": idx_to_a1(row_idx, headers.index(field_name) + 1),
                    "level": "error",
                    "rule ID": f"missing_required_{field_name.lower().replace(' ', '_')}",
                    "rule": f"missing required '{field_name}'",
                    "message": f"add a '{field_name}' term",
                }
            )
        if value and value not in top_terms and value not in valid_labels:
            err_id += 1
            errors.append(
                {
                    "table": table_name,
                    "cell": idx_to_a1(row_idx, headers.index(field_name) + 1),
                    "level": "error",
                    "rule ID": f"invalid_{field_name.lower().replace(' ', '_')}",
                    "rule": f"invalid '{field_name}'",
                    "message": f"replace the '{field_name}' with a term from {source}",
                }
            )
        row_idx += 1
    return errors


def check_restriction_level(table_name, reader, valid_levels):
    """"""
    global err_id
    row_idx = 3
    headers = reader.fieldnames
    errors = []
    for row in reader:
        res_level = row["Restriction Level"]
        if res_level not in valid_levels:
            err_id += 1
            errors.append(
                {
                    "table": table_name,
                    "cell": idx_to_a1(row_idx, headers.index("Restriction Level") + 1),
                    "level": "error",
                    "rule ID": "invalid_restriction_level",
                    "rule": "invalid restriction level",
                    "message": "change the restriction level to one of: "
                    + ", ".join(valid_levels),
                }
            )
        row_idx += 1
    return errors


def idx_to_a1(row, col):
    """Convert a row & column to A1 notation. Adapted from gspread.utils."""
    div = col
    column_label = ""

    while div:
        (div, mod) = divmod(div, 26)
        if mod == 0:
            mod = 26
            div -= 1
        column_label = chr(mod + 64) + column_label

    return f"{column_label}{row}"


def create_message(errors):
    """"""
    table_errors = {}
    for err in errors:
        table_name = err["table"]
        loc = err["cell"]
        row, col = a1_to_idx(loc)
        if table_name in table_errors:
            locs = table_errors[table_name]
        else:
            locs = set()
        locs.add(row)
        table_errors[table_name] = locs

    msg = []
    for table, locs in table_errors.items():
        word = "line"
        if len(locs) > 1:
            word = "lines"
        locs = [str(x) for x in list(locs)]
        rows = ", ".join(locs)
        msg.append(f"- '{table}' {word} {rows}")

    return "\n".join(msg)


def validate_chain(template_dir, labels, genetic_locus_labels, allow_missing):
    """Validate chain.tsv. This checks that:
    - All labels are present in index and end with 'chain'
    - All terms have a parent that is present in this sheet OR is 'protein'
    - All terms with parent 'protein' have a 'Gene' and the 'Gene' is from genetic-locus.tsv

    Return all the labels from this sheet and a set of errors, if any."""
    global err_id
    # Add b2m locus (lives in core but can be used as gene only for b2m)
    table_name = "chain"
    errors = []
    with open(f"{template_dir}/{table_name}.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames
        next(reader)

        # Get the labels for all chains
        # and validate that those labels end with "chain" and are present in index
        if allow_missing:
            chain_labels, label_errors = check_labels(
                table_name,
                reader,
                "index",
                labels,
                regex=r"^.+ chain$",
                missing_level="info",
                allow_missing=allow_missing
            )
        else:
            chain_labels, label_errors = check_labels(
                table_name, reader, "index", labels, regex=r"^.+ chain$"
            )
        errors.extend(label_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate parents
        parent_errors = check_fields(table_name, reader, chain_labels, top_terms=["protein"])
        errors.extend(parent_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)
        row_idx = 3

        # Validate specifics to chain.tsv
        for row in reader:
            # The gene is required when parent == protein
            # The gene should occur in the genetic-locus sheet
            # Don't use check_fields because it is only required when parent == protein
            parent = row["Parent"]
            gene = row["Gene"]
            if not gene or gene.strip == "":
                if parent == "protein":
                    err_id += 1
                    errors.append(
                        {
                            "table": table_name,
                            "cell": idx_to_a1(row_idx, headers.index("Gene") + 1),
                            "level": "error",
                            "rule ID": "missing_chain_gene",
                            "rule": "missing chain gene with 'protein' parent",
                            "message": f"add a 'Gene' from genetic-locus",
                        }
                    )
            elif gene not in genetic_locus_labels:
                err_id += 1
                errors.append(
                    {
                        "table": table_name,
                        "cell": idx_to_a1(row_idx, headers.index("Gene") + 1),
                        "level": "error",
                        "rule ID": "invalid_chain_gene",
                        "rule": "invalid chain gene",
                        "message": f"replace the 'Gene' with a term from genetic-locus",
                    }
                )

            row_idx += 1
    return chain_labels, errors


def validate_chain_sequence(template_dir, chain_labels):
    """Validate chain-sequence.tsv. This checks that:
    - All labels are present in chain.tsv

    Return a set of errors, if any."""
    errors = []
    with open(f"{template_dir}/chain-sequence.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        next(reader)

        # Validate that those labels are present in chain
        # No need to retrieve labels because they're duplicate of those in chain
        _, label_errors = check_labels("chain-sequence", reader, "chain", chain_labels)
        errors.extend(label_errors)
    return errors


def validate_genetic_locus(template_dir, labels, external_labels, allow_missing):
    """Validate genetic-locus.tsv. This checks that:
    - All labels are present in index.tsv and end with 'locus'
    - All terms have parents that are present in this sheet OR is 'MHC locus'
    - Any "In Taxon" value is present in external.tsv

    Return the labels from this sheet and a set of errors, if any."""
    global err_id
    table_name = "genetic-locus"
    errors = []
    with open(f"{template_dir}/{table_name}.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames
        next(reader)

        # Get the labels for all genetic loci
        # and validate that those labels end with "locus" and are present in index
        if allow_missing:
            genetic_locus_labels, label_errors = check_labels(
                table_name,
                reader,
                "index",
                labels,
                regex=r"^.+ locus$",
                missing_level="info",
                allow_missing=allow_missing
            )
        else:
            genetic_locus_labels, label_errors = check_labels(
                table_name, reader, "index", labels, regex=r"^.+ locus$"
            )
        errors.extend(label_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate parents
        parent_errors = check_fields(table_name, reader, genetic_locus_labels, top_terms=["Beta-2-microglobulin locus", "MHC locus"])
        errors.extend(parent_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)
        row_idx = 3

        # Validate taxon
        for row in reader:
            taxon = row["In Taxon"]
            if taxon and taxon.strip() == "":
                taxon = None

            if taxon and taxon not in external_labels:
                err_id += 1
                errors.append(
                    {
                        "table": table_name,
                        "cell": idx_to_a1(row_idx, headers.index("In Taxon") + 1),
                        "level": "error",
                        "rule ID": "invalid_taxon",
                        "rule": "invalid taxon",
                        "message": "add this taxon to 'external' or "
                        "replace it with a taxon from 'external'",
                    }
                )
            row_idx += 1

    return genetic_locus_labels, errors


def validate_halpotype(template_dir, labels, external_labels, allow_missing):
    """Validate haplotype.tsv. This checks that:
    - All labels are present in index.tsv and end with 'haplotype'
    - All terms have a parent that is present in this sheet OR is 'MHC haplotype'
    - An 'In Taxon' value is present when the parent is 'MHC haplotype'
      and this is present in external.tsv

    Return all labels from this sheet and a set of errors, if any."""
    global err_id

    table_name = "haplotype"
    errors = []
    with open(f"{template_dir}/{table_name}.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames
        next(reader)

        # Get the labels for all haplotypes
        # and validate that those labels end with "haplotype" and are present in index
        if allow_missing:
            haplotype_labels, label_errors = check_labels(
                table_name,
                reader,
                "index",
                labels,
                regex=r"^.+ haplotype$",
                missing_level="info",
                allow_missing=allow_missing
            )
        else:
            haplotype_labels, label_errors = check_labels(
                table_name, reader, "index", labels, regex=r"^.+ haplotype$"
            )
        errors.extend(label_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        parent_errors = check_fields(table_name, reader, haplotype_labels, top_terms=["MHC haplotype"])
        errors.extend(parent_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)
        row_idx = 3

        # Validate "In Taxon" values
        for row in reader:
            parent = row["Parent"]
            if parent and parent.strip() == "":
                parent = None
            taxon = row["In Taxon"]
            if taxon and taxon.strip() == "":
                taxon = None

            if parent == "MHC haplotype" and not taxon:
                err_id += 1
                errors.append(
                    {
                        "table": table_name,
                        "cell": idx_to_a1(row_idx, headers.index("In Taxon") + 1),
                        "level": "error",
                        "rule ID": "missing_required_taxon",
                        "rule": "missing required taxon for 'MHC haplotype' parent",
                        "message": "add a taxon from 'external'",
                    }
                )

            if taxon and taxon not in external_labels:
                err_id += 1
                errors.append(
                    {
                        "table": table_name,
                        "cell": idx_to_a1(row_idx, headers.index("In Taxon") + 1),
                        "level": "error",
                        "rule ID": "invalid_taxon",
                        "rule": "invalid taxon",
                        "message": "add this taxon to 'external' or "
                        "replace it with a taxon from 'external'",
                    }
                )
            row_idx += 1
    return haplotype_labels, errors


def validate_haplotype_molecule(
    template_dir, labels, molecule_labels, haplotype_labels, external_labels, allow_missing
):
    """Validate haplotype-molecule.tsv. This checks that:
    - All labels are present in index.tsv and end with 'with X haplotype' or 'with haplotype'
    - All terms have a parent that is present in molecule.tsv
    - The 'Restriction Level' is 'haplotype'
    - All terms have an "In Taxon" value that is present in external.tsv
    - All terms have a "With Haplotype" value that is present in haplotype.tsv

    Return a set of errors, if any."""
    global err_id
    table_name = "haplotype-molecule"
    errors = []
    with open(f"{template_dir}/{table_name}.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames
        next(reader)

        # Get any labels not defined in index
        if allow_missing:
            _, label_errors = check_labels(
                table_name,
                reader,
                "index",
                labels,
                regex=r"^.+ (with haplotype|with [^ ]+ haplotype)$",
                missing_level="info",
                allow_missing=allow_missing
            )
        else:
            _, label_errors = check_labels(
                table_name,
                reader,
                "index",
                labels,
                regex=r"^.+ (with haplotype|with [^ ]+ haplotype)$",
            )
        errors.extend(label_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate parents (from molecule table)
        parent_errors = check_fields(
            table_name, reader, molecule_labels, top_terms=["MHC protein complex"], source="molecule",
        )
        errors.extend(parent_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate with haplotype
        with_haplotype_errors = check_fields(
            table_name,
            reader,
            haplotype_labels,
            top_terms=["MHC haplotype"],
            field_name="With Haplotype",
            source="haplotype",
        )
        errors.extend(with_haplotype_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate restriction levels
        res_level_errors = check_restriction_level(table_name, reader, ["haplotype"])
        errors.extend(res_level_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)
        row_idx = 3

        # Validate taxon
        for row in reader:
            taxon = row["In Taxon"]
            if taxon.strip == "":
                taxon = None

            if not taxon:
                err_id += 1
                errors.append(
                    {
                        "table": table_name,
                        "cell": idx_to_a1(row_idx, headers.index("In Taxon") + 1),
                        "level": "error",
                        "rule ID": "missing_required_taxon",
                        "rule": "missing required taxon",
                        "message": "add a taxon from 'external'",
                    }
                )

            elif taxon not in external_labels:
                err_id += 1
                errors.append(
                    {
                        "table": table_name,
                        "cell": idx_to_a1(row_idx, headers.index("In Taxon") + 1),
                        "level": "error",
                        "rule ID": "invalid_taxon",
                        "rule": "invalid taxon",
                        "message": "add this taxon to 'external' or "
                        "replace it with a taxon from 'external'",
                    }
                )
            row_idx += 1
    return errors


def validate_iedb_labels(iedb_path, labels, allow_missing):
    """Validate iedb.tsv. This checks that all labels are present in index.tsv. Return a set of
    errors, if any."""
    errors = []
    with open(iedb_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        next(reader)
        if not allow_missing:
            _, errors = check_labels("iedb", reader, "index", labels)
    return errors


def validate_molecule(
    template_dir,
    labels,
    chain_labels,
    external_labels,
    haplotype_labels,
    serotype_labels,
    allow_missing,
):
    """Validate molecule.tsv. This checks that:
    - All labels are present in index.tsv and end with 'protein complex'
    - All terms have parents that are present in this sheet OR is 'MHC protein complex'
    - The 'Restriction Level' is one of: class, complete molecule, locus, partial molecule
    - An 'In Taxon' value is present when the parent is NOT 'MHC protein complex'
      and this value is present in external.tsv
    - The 'Alpha Chain' value is present in chain.tsv
    - The 'Beta Chain' value is present in chain.tsv OR is 'Beta-2-microglobulin'
    - The 'With Haplotype' value is present in haplotype.tsv
    - The 'With Serotype' value is present in serotype.tsv

    Return all labels from this sheet and a set of errors, if any."""
    global err_id
    table_name = "molecule"
    errors = []
    with open(f"{template_dir}/molecule.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames
        next(reader)

        # Get the labels for all genetic loci
        # and validate that those labels end with "locus" and are present in index
        if allow_missing:
            molecule_labels, label_errors = check_labels(
                "molecule",
                reader,
                "index",
                labels,
                regex=r"^.+ protein complex$",
                missing_level="info",
                allow_missing=allow_missing
            )
        else:
            molecule_labels, label_errors = check_labels(
                "molecule", reader, "index", labels, regex=r"^.+ protein complex$"
            )
        errors.extend(label_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate parents
        parent_errors = check_fields(
            table_name, reader, molecule_labels, top_terms=["MHC protein complex"]
        )
        errors.extend(parent_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Alpha chain must be in chain labels
        alpha_errors = check_fields(
            table_name,
            reader,
            chain_labels,
            field_name="Alpha Chain",
            source="chain",
            required=False,
        )
        errors.extend(alpha_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Beta chain must be in chain labels or Beta-2-microglobulin
        beta_errors = check_fields(
            table_name,
            reader,
            chain_labels,
            field_name="Beta Chain",
            top_terms=["Beta-2-microglobulin"],
            source="chain",
            required=False,
        )
        errors.extend(beta_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # With haplotype is in haplotype
        haplotype_errors = check_fields(
            table_name,
            reader,
            haplotype_labels,
            field_name="With Haplotype",
            source="haplotype",
            required=False,
        )
        errors.extend(haplotype_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # With serotype is in serotype
        serotype_errors = check_fields(
            table_name,
            reader,
            serotype_labels,
            field_name="With Serotype",
            source="serotype",
            required=False,
        )
        errors.extend(serotype_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate restriction levels
        res_level_errors = check_restriction_level(
            table_name, reader, ["class", "locus", "complete molecule", "partial molecule",],
        )
        errors.extend(res_level_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)
        row_idx = 3

        # Restriction level must be one of: locus, complete molecule, partial molecule
        # Taxon required when parent NOT MHC protein complex
        for row in reader:
            parent = row["Parent"]
            taxon = row["In Taxon"]
            if taxon.strip == "":
                taxon = None

            if parent != "MHC protein complex" and not taxon:
                err_id += 1
                errors.append(
                    {
                        "table": table_name,
                        "cell": idx_to_a1(row_idx, headers.index("In Taxon") + 1),
                        "level": "error",
                        "rule ID": "missing_required_taxon",
                        "rule": "missing required taxon for parent other than "
                        "'MHC protein complex'",
                        "message": "add a taxon from 'external'",
                    }
                )

            elif taxon and taxon not in external_labels:
                err_id += 1
                errors.append(
                    {
                        "table": table_name,
                        "cell": idx_to_a1(row_idx, headers.index("In Taxon") + 1),
                        "level": "error",
                        "rule ID": "invalid_taxon",
                        "rule": "invalid taxon",
                        "message": "add this taxon to 'external' or "
                        "replace it with a taxon from 'external'",
                    }
                )
            row_idx += 1

    return molecule_labels, errors


def validate_mutant_molecule(template_dir, labels, external_labels, molecule_labels, allow_missing):
    """Validate mutant-molecule.tsv. This checks that:
    - All labels are present in index.tsv and end with 'protein complex'
    - All terms have parents that are present in this sheet OR is 'mutant MHC protein complex'
    - The 'Restriction Level' is one of: class, complete molecule, partial molecule
    - The 'In Taxon' value is in external.tsv
    - The 'Mutant Of' value is in molecule.tsv

    Return a set of errors, if any."""
    global err_id
    table_name = "mutant-molecule"
    errors = []
    with open(f"{template_dir}/{table_name}.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames
        next(reader)

        # Get any labels not defined in index
        if allow_missing:
            mutant_molecule_labels, label_errors = check_labels(
                table_name,
                reader,
                "index",
                labels,
                regex=r"^.+ protein complex$",
                missing_level="info",
                allow_missing=allow_missing
            )
        else:
            mutant_molecule_labels, label_errors = check_labels(
                table_name, reader, "index", labels, regex=r"^.+ protein complex$"
            )
        errors.extend(label_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate parents
        parent_errors = check_fields(
            table_name, reader, mutant_molecule_labels, top_terms=["mutant MHC protein complex"],
        )
        errors.extend(parent_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate "Mutant Of" fields
        molecule_errors = check_fields(
            table_name,
            reader,
            molecule_labels,
            field_name="Mutant Of",
            source="molecule",
            required=False,
        )
        errors.extend(molecule_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate restriction levels
        res_level_errors = check_restriction_level(
            table_name, reader, ["class", "complete molecule", "partial molecule"]
        )
        errors.extend(res_level_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)
        row_idx = 3

        # Taxon required when parent NOT MHC protein complex
        for row in reader:
            taxon = row["In Taxon"]
            if taxon.strip == "":
                taxon = None

            if taxon and taxon not in external_labels:
                err_id += 1
                errors.append(
                    {
                        "table": table_name,
                        "cell": idx_to_a1(row_idx, headers.index("In Taxon") + 1),
                        "level": "error",
                        "rule ID": "invalid_taxon",
                        "rule": "invalid taxon",
                        "message": "add this taxon to 'external' or "
                        "replace it with a taxon from 'external'",
                    }
                )
            row_idx += 1
    return errors


def validate_serotype(template_dir, labels, external_labels, allow_missing):
    """Validate serotype.tsv. This checks that:
    - All labels are present in index.tsv and end with 'serotype'
    - All terms have parents that are present in this sheet OR is 'MHC serotype'
    - The 'In Taxon' value is present in external.tsv

    Return all labels from this sheet and a set of errors, if any.
    - """
    global err_id
    table_name = "serotype"
    errors = []
    with open(f"{template_dir}/{table_name}.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames
        next(reader)

        # Get the labels for all serotypes
        # and validate that those labels end with "serotype" and are present in index
        if allow_missing:
            serotype_labels, label_errors = check_labels(
                table_name,
                reader,
                "index",
                labels,
                regex=r"^.+ serotype$",
                missing_level="info",
                allow_missing=allow_missing
            )
        else:
            serotype_labels, label_errors = check_labels(
                table_name, reader, "index", labels, regex=r"^.+ serotype$"
            )
        errors.extend(label_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate parents
        parent_errors = check_fields(table_name, reader, serotype_labels, top_terms=["MHC serotype"])
        errors.extend(parent_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)
        row_idx = 3

        # Validate taxon in 'external'
        for row in reader:
            taxon = row["In Taxon"]
            if taxon and taxon.strip() == "":
                taxon = None

            if taxon and taxon not in external_labels:
                err_id += 1
                errors.append(
                    {
                        "table": table_name,
                        "cell": idx_to_a1(row_idx, headers.index("In Taxon") + 1),
                        "level": "error",
                        "rule ID": "invalid_taxon",
                        "rule": "invalid taxon",
                        "message": "add this taxon to 'external' or "
                        "replace it with a taxon from 'external'",
                    }
                )
            row_idx += 1
    return serotype_labels, errors


def validate_serotype_molecule(
    template_dir, labels, external_labels, molecule_labels, serotype_labels, allow_missing
):
    """Validate serotype-molecule.tsv. This checks that:
    - All labels are present in index.tsv and end with 'with X serotype' or 'with serotype'
    - All terms have parents that are present in molecule.tsv
    - The 'Restriction Level' value is 'serotype'
    - The 'In Taxon' value is present in external.tsv
    - All terms have 'With Serotype' values that are present in serotype.tsv

    Return a set of errors, if any."""
    global err_id
    table_name = "serotype-molecule"
    errors = []
    with open(f"{template_dir}/{table_name}.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames
        next(reader)

        if allow_missing:
            _, label_errors = check_labels(
                table_name,
                reader,
                "index",
                labels,
                regex=r"^.+ (with serotype|with [^ ]+ serotype)$",
                missing_level="info",
                allow_missing=allow_missing
            )
        else:
            _, label_errors = check_labels(
                table_name,
                reader,
                "index",
                labels,
                regex=r"^.+ (with serotype|with [^ ]+ serotype)$",
            )
        errors.extend(label_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate parents (from molecule)
        parent_errors = check_fields(
            table_name, reader, molecule_labels, top_terms=["MHC protein complex"], source="molecule",
        )
        errors.extend(parent_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        with_serotype_errors = check_fields(
            table_name, reader, serotype_labels, field_name="With Serotype", source="serotype",
        )
        errors.extend(with_serotype_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)

        # Validate restriction levels
        res_level_errors = check_restriction_level(table_name, reader, ["serotype"])
        errors.extend(res_level_errors)

        # Reset file to beginning
        f.seek(0)
        next(reader)
        next(reader)
        row_idx = 3

        # Validate taxon
        for row in reader:
            taxon = row["In Taxon"]
            if taxon.strip == "":
                taxon = None

            if taxon and taxon not in external_labels:
                err_id += 1
                errors.append(
                    {
                        "table": table_name,
                        "cell": idx_to_a1(row_idx, headers.index("In Taxon") + 1),
                        "level": "error",
                        "rule ID": "invalid_taxon",
                        "rule": "invalid taxon",
                        "message": "add this taxon to 'external' or "
                        "replace it with a taxon from 'external'",
                    }
                )
            row_idx += 1
    return errors


def main():
    p = ArgumentParser()
    p.add_argument("index")
    p.add_argument("iedb")
    p.add_argument("template_dir")
    p.add_argument("err_output")
    p.add_argument("-a", "--allow-missing", action="store_true")
    args = p.parse_args()

    template_dir = args.template_dir
    allow_missing = args.allow_missing

    # Get MRO labels defined in the index
    labels = []
    with open(args.index, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        next(reader)
        for row in reader:
            l = row["Label"]
            if l and l.strip() != "":
                labels.append(l)

    # Get imported term labels
    ext_labels = []
    with open(f"{template_dir}/external.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        next(reader)
        for row in reader:
            l = row["Label"]
            if l and l.strip() != "":
                ext_labels.append(l)

    errors = []

    # Validate the IEDB table
    iedb_errors = validate_iedb_labels(args.iedb, labels, allow_missing)
    errors.extend(iedb_errors)

    # Validate genetic-locus
    genetic_locus_labels, genetic_locus_errors = validate_genetic_locus(
        template_dir, labels, ext_labels, allow_missing
    )

    # Validate chain
    chain_labels, chain_errors = validate_chain(
        template_dir, labels, genetic_locus_labels, allow_missing
    )

    # Validate chain-sequence
    chain_sequence_errors = validate_chain_sequence(template_dir, chain_labels)

    # Vaildate haplotype
    haplotype_labels, haplotype_errors = validate_halpotype(
        template_dir, labels, ext_labels, allow_missing
    )

    # Validate serotype
    serotype_labels, serotype_errors = validate_serotype(
        template_dir, labels, ext_labels, allow_missing
    )

    # Validate molecule
    molecule_labels, molecule_errors = validate_molecule(
        template_dir,
        labels,
        chain_labels,
        ext_labels,
        haplotype_labels,
        serotype_labels,
        allow_missing,
    )

    # Validate mutant-molecule
    mutant_molecule_errors = validate_mutant_molecule(
        template_dir, labels, ext_labels, molecule_labels, allow_missing
    )

    # Validate haplotype-moleculee
    haplotype_molecule_errors = validate_haplotype_molecule(
        template_dir, labels, molecule_labels, haplotype_labels, ext_labels, allow_missing
    )

    # Validate serotype-molecule
    serotype_molecule_errors = validate_serotype_molecule(
        template_dir, labels, ext_labels, molecule_labels, serotype_labels, allow_missing
    )

    # Add errors in table order
    errors.extend(chain_errors)
    errors.extend(chain_sequence_errors)
    errors.extend(genetic_locus_errors)
    errors.extend(haplotype_errors)
    errors.extend(haplotype_molecule_errors)
    errors.extend(molecule_errors)
    errors.extend(mutant_molecule_errors)
    errors.extend(serotype_errors)
    errors.extend(serotype_molecule_errors)

    with open(args.err_output, "w") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "table",
                "cell",
                "level",
                "rule ID",
                "rule",
                "message",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(errors)

    if errors:
        has_error_level = [x for x in errors if x["level"] == "error"]
        msg = create_message(errors)
        print(
            "\n-------------------------------------------------------\n"
            f"ERROR: Validation completed with {err_id} message(s) at:"
            f"\n{msg}\n"
            "---------------------------------------------------------\n"
        )
    else:
        print("\nValidation passed!\n")


if __name__ == "__main__":
    main()
