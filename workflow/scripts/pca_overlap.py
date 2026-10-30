"""List reference PCA variant IDs (.eigenvec.allele) present in a .pvar file.

Fails if fewer than --min-overlap of the reference variants are found.
"""

import argparse
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pvar", required=True)
    p.add_argument("--eigenvec", required=True)
    p.add_argument("--min-overlap", type=float, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--summary", required=True)
    args = p.parse_args()

    with open(args.eigenvec) as f:
        next(f)
        ref = {line.split("\t", 2)[1] for line in f}

    found = set()
    with open(args.pvar) as f:
        for line in f:
            if line.startswith("#"):
                continue
            vid = line.split("\t", 3)[2]
            if vid in ref:
                found.add(vid)

    with open(args.out, "w") as f:
        f.writelines(f"{v}\n" for v in sorted(found))

    frac = len(found) / len(ref)
    with open(args.summary, "w") as f:
        f.write(f"pca_ref_variants\t{len(ref)}\n")
        f.write(f"pca_ref_variants_found\t{len(found)}\n")
        f.write(f"pca_ref_fraction_found\t{frac:.4f}\n")
    print(f"{len(found)} of {len(ref)} reference PCA variants found ({frac:.2%})")

    if frac < args.min_overlap:
        sys.exit(f"fraction of reference PCA variants found below {args.min_overlap}")


if __name__ == "__main__":
    main()
