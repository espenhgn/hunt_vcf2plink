"""Lift .bim variant positions with UCSC liftOver and flag sites in inverted regions.

Implements the triple-liftOver method (Sheng et al. 2022, doi:10.1101/2022.02.19.481174):
each variant is lifted together with the preceding and succeeding base. The
variant is "inverted" if the preceding base maps right after it in the target
build, or the succeeding base right before it; SNV alleles there must be strand
flipped. Variants that do not lift, or lift to another chromosome, are left out.

Output (tab-separated, with header): snpid, chr, target position, base position,
category ("lifted" or "inverted").
"""

import argparse
import subprocess
import sys
from pathlib import Path

# PLINK numeric/alias chromosome codes -> UCSC names
CHR_ALIASES = {"23": "X", "25": "X", "XY": "X", "PAR1": "X", "PAR2": "X",
               "24": "Y", "26": "M", "MT": "M"}


def ucsc_chr(c):
    c = c.removeprefix("chr")
    return CHR_ALIASES.get(c, c)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bim", required=True)
    p.add_argument("--liftover", required=True, help="UCSC liftOver binary")
    p.add_argument("--chain", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    bed, lifted_bed, unmapped = (work / n for n in ("in.bed", "lifted.bed", "unmapped.bed"))

    # BED (0-based, half-open) entries: preceding (a), variant (b), succeeding (c) base
    variants = []  # (id, chr, pos)
    with open(args.bim) as f, open(bed, "w") as out:
        for i, line in enumerate(f):
            c, vid, _, pos = line.split()[:4]
            c, pos = ucsc_chr(c), int(pos)
            variants.append((vid, c, pos))
            out.write(f"chr{c}\t{pos - 1}\t{pos}\t{i}_b\n")
            if pos > 1:
                out.write(f"chr{c}\t{pos - 2}\t{pos - 1}\t{i}_a\n")
            out.write(f"chr{c}\t{pos}\t{pos + 1}\t{i}_c\n")

    res = subprocess.run([args.liftover, bed, args.chain, lifted_bed, unmapped],
                         capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"liftOver failed ({res.returncode}):\n{res.stdout}{res.stderr}")

    # lifted[tag][i] = (chr, 1-based position)
    lifted = {"a": {}, "b": {}, "c": {}}
    with open(lifted_bed) as f:
        for line in f:
            c, _, end, name = line.split("\t")[:4]
            i, tag = name.rstrip("\n").split("_")
            lifted[tag][int(i)] = (c.removeprefix("chr"), int(end))

    n_out = n_inverted = 0
    with open(args.out, "w") as out:
        out.write("snpid\tchr\tpos_target\tpos_base\tcategory\n")
        for i, (vid, c, pos) in enumerate(variants):
            if i not in lifted["b"] or lifted["b"][i][0] != c:
                continue
            new_pos = lifted["b"][i][1]
            prev, succ = lifted["a"].get(i), lifted["c"].get(i)
            inverted = (prev is not None and prev[0] == c and prev[1] - new_pos == 1) or \
                       (succ is not None and succ[0] == c and succ[1] - new_pos == -1)
            n_inverted += inverted
            n_out += 1
            out.write(f"{vid}\t{c}\t{new_pos}\t{pos}\t{'inverted' if inverted else 'lifted'}\n")

    print(f"{len(variants)} variants in, {n_out} lifted to the same chromosome, "
          f"{n_inverted} inverted")


if __name__ == "__main__":
    main()
