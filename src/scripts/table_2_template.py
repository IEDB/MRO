import csv
import os
import json
import sqlite3

from argparse import ArgumentParser


def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d


def main():
    parser = ArgumentParser()
    parser.add_argument("db")
    parser.add_argument("template_dir")
    parser.add_argument("name")
    args = parser.parse_args()

    template = os.path.join(args.template_dir, args.name + ".tsv")

    # Get the current order of the headers
    with open(template, "r") as f:
        headers = f.readline().strip().split("\t")
    print(headers)

    with sqlite3.connect(args.db) as conn:
        conn.row_factory = dict_factory
        cur = conn.cursor()

        # Get the ROBOT template strings (we always use an underscore instead of dash in db)
        results = cur.execute(
            'SELECT column, template FROM "column" WHERE "table" = ?',
            (args.name.replace("-", "_"),)
        ).fetchall()
        robot_templates = {res["column"]: res["template"] for res in results}

        # Get the rows of the template, using the meta columns for invalid values
        rows = [robot_templates]
        results = cur.execute(f'SELECT * FROM "{args.name.replace("-", "_")}_view"').fetchall()
        for res in results:
            row = {}
            for col, val in res.items():
                if col == "row_number" or col.endswith("_meta"):
                    continue
                meta_val = res[col + "_meta"]
                if not meta_val:
                    row[col] = val
                    continue
                meta_val = json.loads(meta_val)
                row[col] = meta_val["value"]
            rows.append(row)

    with open(template, "w") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
