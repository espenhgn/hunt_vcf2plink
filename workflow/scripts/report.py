"""Collect per-step variant counts into one TSV (one row per chromosome, plus total).

Sample and PCA summaries are written as leading '#' comment lines.
"""

import argparse
import re

B37_PATTERNS = {
    "vcf_variants": re.compile(r"--vcf: (\d+) variants scanned"),
    "b37_after_qc": re.compile(r"^(\d+) variants remaining after main filters", re.M),
}


def read_kv(path):
    with open(path) as f:
        return dict(line.rstrip("\n").split("\t", 1) for line in f if line.strip())


def read_b37_log(path):
    text = open(path).read()
    out = {}
    for key, pat in B37_PATTERNS.items():
        m = pat.findall(text)
        out[key] = m[-1] if m else "NA"
    return out


def count_variants(pvar):
    with open(pvar) as f:
        return sum(1 for line in f if not line.startswith("#"))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sex", required=True)
    p.add_argument("--liftover", nargs="+", required=True)
    p.add_argument("--align", nargs="+", required=True)
    p.add_argument("--b37-logs", nargs="+", required=True)
    p.add_argument("--merged-pvar", required=True)
    p.add_argument("--pca", required=True)
    args = p.parse_args()

    for path in (args.sex, args.pca):
        for k, v in read_kv(path).items():
            print(f"# {k}\t{v}")
    print(f"# merged_aligned_variants\t{count_variants(args.merged_pvar)}")

    rows = []
    for b37_log, lift, align in zip(args.b37_logs, args.liftover, args.align):
        row = read_kv(lift)
        row.update(read_b37_log(b37_log))
        row.update({k: v for k, v in read_kv(align).items() if k != "chr"})
        rows.append(row)

    cols = ["chr", "vcf_variants", "b37_after_qc"]
    cols += [k for k in rows[0] if k not in cols]
    total = {"chr": "total"}
    for c in cols[1:]:
        vals = [r.get(c, "NA") for r in rows]
        total[c] = str(sum(int(v) for v in vals)) if all(v.isdigit() for v in vals) else "NA"

    print("\t".join(cols))
    for r in rows + [total]:
        print("\t".join(r.get(c, "NA") for c in cols))


if __name__ == "__main__":
    main()
