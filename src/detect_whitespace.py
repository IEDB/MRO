import csv
import re
import sys

from argparse import ArgumentParser


def main():
    p = ArgumentParser()
    p.add_argument("paths", nargs="+")
    p.add_argument("output")
    args = p.parse_args()

    errs = []
    err_tables = 0
    for path in args.paths:
        has_err = False
        with open(path, "r") as f:
            if path.endswith(".csv"):
                reader = csv.reader(f)
            else:
                reader = csv.reader(f, delimiter="\t")

            row_num = 1
            for row in reader:
                col_num = 1
                for cell in row:
                    if re.findall(r"^\s+", cell) or re.findall(r"\s+$", cell):
                        has_err = True
                        errs.append(
                            {
                                "table": path,
                                "column": col_num,
                                "row": row_num,
                                "content": cell,
                            }
                        )
                    col_num += 1
                row_num += 1
        if has_err:
            err_tables += 1

    if errs:
        with open(args.output, "w") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["table", "column", "row", "content"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(errs)
        print(
            f"ERROR: {err_tables} table(s) have one or more cells with leading or trailing whitespace!"
            f"\n       See {args.output} for details."
        )
        sys.exit(1)


if __name__ == '__main__':
    main()
