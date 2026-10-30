"""Write a plink2 --update-sex file (#IID, SEX) from the covariates file.

Only sex codes 1 (male) and 2 (female) are written. Fails if fewer than
--min-match of the VCF samples get a sex, which catches ID column mismatches.
"""

import argparse
import csv
import gzip
import sys


def open_text(path):
    with open(path, "rb") as f:
        gzipped = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gzipped else open(path)


def vcf_samples(path):
    with open_text(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    sys.exit(f"no #CHROM header line in {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--covariates", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--id-col", required=True)
    p.add_argument("--sex-col", required=True)
    p.add_argument("--min-match", type=float, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--summary", required=True)
    args = p.parse_args()

    sex = {}
    with open_text(args.covariates) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for col in (args.id_col, args.sex_col):
            if col not in reader.fieldnames:
                sys.exit(f"column {col!r} not in {args.covariates}: {reader.fieldnames}")
        for row in reader:
            if row[args.sex_col] in ("1", "2"):
                sex[row[args.id_col]] = row[args.sex_col]

    samples = vcf_samples(args.vcf)
    matched = [s for s in samples if s in sex]
    n_male = sum(sex[s] == "1" for s in matched)
    frac = len(matched) / len(samples) if samples else 0.0

    with open(args.out, "w") as f:
        f.write("#IID\tSEX\n")
        for s in matched:
            f.write(f"{s}\t{sex[s]}\n")

    stats = {
        "covariate_rows_with_sex": len(sex),
        "vcf_samples": len(samples),
        "vcf_samples_with_sex": len(matched),
        "males": n_male,
        "females": len(matched) - n_male,
        "fraction_with_sex": f"{frac:.4f}",
    }
    with open(args.summary, "w") as f:
        for k, v in stats.items():
            f.write(f"{k}\t{v}\n")
            print(f"{k}\t{v}")

    if frac < args.min_match:
        sys.exit(
            f"only {frac:.2%} of VCF samples have sex in {args.covariates} "
            f"(min {args.min_match:.0%}); check COVARIATES_ID_COL"
        )


if __name__ == "__main__":
    main()
