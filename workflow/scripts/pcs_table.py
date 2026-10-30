"""Write IID and PC1..PCn from a plink2 --score .sscore file as TSV to stdout.

usage: pcs_table.py SSCORE N_PC
"""

import csv
import sys


def main():
    sscore, n_pc = sys.argv[1], int(sys.argv[2])
    with open(sscore) as f:
        reader = csv.reader(f, delimiter="\t")
        header = [h.lstrip("#") for h in next(reader)]
        iid = header.index("IID")
        # with variance-standardize, score columns are named PC1_AVG, ...
        pcs = [header.index(f"PC{i}_AVG") for i in range(1, n_pc + 1)]
        out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        out.writerow(["IID"] + [f"PC{i}" for i in range(1, n_pc + 1)])
        for row in reader:
            out.writerow([row[iid]] + [row[i] for i in pcs])


if __name__ == "__main__":
    main()
