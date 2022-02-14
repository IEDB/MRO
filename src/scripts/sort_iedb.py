import csv

from argparse import ArgumentParser


def main():
    parser = ArgumentParser()
    parser.add_argument("iedb")
    args = parser.parse_args()

    with open(args.iedb, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        robot_header = next(reader)
        rows_sorted = sorted(reader, key=lambda d: int(d["IEDB ID"]))
    with open(args.iedb, "w") as f:
        writer = csv.DictWriter(
            f, delimiter="\t", fieldnames=rows_sorted[0].keys(), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(robot_header)
        writer.writerows(rows_sorted)


if __name__ == '__main__':
    main()
