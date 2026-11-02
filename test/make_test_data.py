"""Generate a small synthetic test data set mimicking the HUNT VCF files (GRCh37).

- chr22: the reference PCA SNPs (GRCh38) lifted to GRCh37, with genotypes from two
  pseudo-populations, a fraction with swapped or strand-flipped alleles, plus indels.
- chrX: random SNPs in PAR1, non-PAR and PAR2 (GRCh37); males haploid outside PAR.
- A GRCh38 FASTA (chr22, chrX) that is N except at the test SNP positions.
- A covariates file with PID and Sex columns.
- A PCA basis (.acount, .eigenvec.allele, .eigenval) from a simulated GRCh38
  reference panel of the same two pseudo-populations, made with plink2 --pca.

No real data is used. Outputs are committed to test/data; rerun only to change them.

The chr22 sites are taken from the HGDP+1KG PCA reference (.acount), which is not
part of this repository.

usage: python make_test_data.py --liftover liftOver --plink2 plink2
           --hg38-to-hg19 hg38ToHg19.over.chain.gz --hg19-to-hg38 hg19ToHg38.over.chain.gz
           --acount hgdp1kgp_chr_autosomes.pca_unrelated.acount --out data
"""

import argparse
import gzip
import random
import shutil
import subprocess
import tempfile
from pathlib import Path

N_SAMPLES = 100
N_PANEL = 200
N_PC = 10
SEED = 20261030
CPL = {"A": "T", "T": "A", "C": "G", "G": "C"}
# GRCh37 contig lengths, as in the HUNT VCF header
B37_LEN = {"22": 51304566, "X": 155270560}
# GRCh38 contig lengths
B38_LEN = {"22": 50818468, "X": 156040895}


def liftover(liftover_bin, chain, chrom, positions):
    """Return {pos: lifted_pos} for positions that lift to the same chromosome."""
    with tempfile.TemporaryDirectory() as d:
        bed, out, unmapped = (Path(d) / n for n in ("in.bed", "out.bed", "unmapped.bed"))
        bed.write_text("".join(f"chr{chrom}\t{p - 1}\t{p}\t{p}\n" for p in positions))
        subprocess.run([liftover_bin, bed, chain, out, unmapped], check=True, capture_output=True)
        lifted = {}
        for line in out.read_text().splitlines():
            c, _, end, name = line.split("\t")[:4]
            if c == f"chr{chrom}":
                lifted[int(name)] = int(end)
        return lifted


def random_snp_alleles(rng):
    while True:
        ref, alt = rng.sample("ACGT", 2)
        if CPL[ref] != alt:
            return ref, alt


def genotypes(rng, freqs, haploid):
    """VCF GT:DS fields; freqs per sample (pseudo-population), haploid flags per sample."""
    fields = []
    for p, hap in zip(freqs, haploid):
        if hap:
            g = int(rng.random() < p)
            ds = min(max(g + rng.gauss(0, 0.03), 0), 1)
            fields.append(f"{g}:{ds:.3f}")
        else:
            a, b = int(rng.random() < p), int(rng.random() < p)
            ds = min(max(a + b + rng.gauss(0, 0.03), 0), 2)
            fields.append(f"{a}|{b}:{ds:.3f}")
    return fields


def vcf_header(chrom, samples):
    return (
        "##fileformat=VCFv4.2\n"
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        f"##contig=<ID={chrom},length={B37_LEN[chrom]}>\n"
        '##INFO=<ID=INFO,Number=1,Type=Float,Description="IMPUTE2 info score">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Genotype dosage">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
    )


def write_vcf(path, chrom, samples, records):
    with gzip.open(path, "wt") as f:
        f.write(vcf_header(chrom, samples))
        for pos, ref, alt, info, gts in sorted(records):
            f.write(f"{chrom}\t{pos}\t{chrom}:{pos}_{ref}/{alt}\t{ref}\t{alt}\t.\tPASS\t"
                    f"INFO={info:.3f}\tGT:DS\t" + "\t".join(gts) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--liftover", required=True)
    p.add_argument("--plink2", required=True)
    p.add_argument("--hg38-to-hg19", required=True)
    p.add_argument("--hg19-to-hg38", required=True)
    p.add_argument("--acount", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    rng = random.Random(SEED)
    out = Path(args.out)
    for d in ("vcf", "ref", "covariates", "pca"):
        (out / d).mkdir(parents=True, exist_ok=True)

    samples = [f"S{i:03d}" for i in range(1, N_SAMPLES + 1)]
    sex = {s: 1 + i % 2 for i, s in enumerate(samples)}  # alternating male/female
    pop2 = {s: i >= N_SAMPLES // 2 for i, s in enumerate(samples)}
    fasta = {c: {} for c in B38_LEN}  # chrom -> {b38 pos: base}

    def info(rng):
        return rng.uniform(0.1, 0.3) if rng.random() < 0.05 else rng.uniform(0.6, 1.0)

    # chr22: reference PCA SNPs
    ref_snps = []
    with open(args.acount) as f:
        next(f)
        for line in f:
            c, _, ref, alt, alt_ct, obs_ct = line.split()
            vid = line.split()[1]
            if c == "22":
                ref_snps.append((int(vid.split(":")[1]), ref, alt, int(alt_ct) / int(obs_ct)))
    to_b37 = liftover(args.liftover, args.hg38_to_hg19, "22", [s[0] for s in ref_snps])
    records = []
    sites = []  # (b38 pos, ref, alt, pop1 freq, pop2 freq), for the PCA reference panel
    for pos38, ref, alt, freq in ref_snps:
        if pos38 not in to_b37:
            continue
        fasta["22"][pos38] = ref
        p2 = min(max(freq + rng.gauss(0, 0.2), 0.02), 0.98)
        sites.append((pos38, ref, alt, freq, p2))
        freqs = [p2 if pop2[s] else freq for s in samples]
        a, b = ref, alt
        u = rng.random()
        if u < 0.10:  # REF/ALT swapped relative to the reference sequence
            a, b, freqs = b, a, [1 - q for q in freqs]
        elif u < 0.15:  # opposite strand
            a, b = CPL[a], CPL[b]
        records.append((to_b37[pos38], a, b, info(rng), genotypes(rng, freqs, [False] * N_SAMPLES)))
    used = {r[0] for r in records}
    for _ in range(20):  # indels, kept only in the all-variant GRCh38 set
        pos = rng.randrange(20_000_000, 50_000_000)
        if pos not in used:
            ref = rng.choice("ACGT")
            records.append((pos, ref, ref + rng.choice("ACGT"), 0.9,
                            genotypes(rng, [0.3] * N_SAMPLES, [False] * N_SAMPLES)))
    write_vcf(out / "vcf" / "CHR22_TEST.vcf.gz", "22", samples, records)

    # chrX: PAR1 (b37 <= 2699520), non-PAR, PAR2 (b37 >= 154931044)
    regions = [(200_000, 2_600_000, 100, True), (3_000_000, 150_000_000, 200, False),
               (154_940_000, 155_200_000, 20, True)]
    candidates = []
    for lo, hi, n, par in regions:
        candidates += [(pos, par) for pos in rng.sample(range(lo, hi), n)]
    to_b38 = liftover(args.liftover, args.hg19_to_hg38, "X", [c[0] for c in candidates])
    records = []
    for pos37, par in candidates:
        if pos37 not in to_b38:
            continue
        ref, alt = random_snp_alleles(rng)
        fasta["X"][to_b38[pos37]] = ref
        freq = rng.uniform(0.05, 0.5)
        haploid = [sex[s] == 1 and not par for s in samples]
        records.append((pos37, ref, alt, info(rng), genotypes(rng, [freq] * N_SAMPLES, haploid)))
    write_vcf(out / "vcf" / "CHRX_TEST.vcf.gz", "X", samples, records)

    # GRCh38 FASTA, N except at test SNPs; compressed with zstd by the caller
    with open(out / "ref" / "GRCh38_test.fa", "w") as f:
        for c, length in B38_LEN.items():
            seq = bytearray(b"N" * length)
            for pos, base in fasta[c].items():
                seq[pos - 1] = ord(base)
            f.write(f">chr{c}\n")
            for i in range(0, length, 60):
                f.write(seq[i:i + 60].decode() + "\n")

    # covariates: last sample has no sex; one extra row not in the VCF
    with gzip.open(out / "covariates" / "covariates_TEST.txt.gz", "wt") as f:
        f.write("PID\tcurrent\tImputed20220624\tBatchDetailed\tSex\tANCESTRY\n")
        for s in samples:
            f.write(f"{s}\t1\t1\tbatch1\t{'NA' if s == samples[-1] else sex[s]}\tEUR\n")
        f.write("S999\t1\t1\tbatch1\t2\tEUR\n")

    make_pca_reference(args.plink2, sites, out / "pca" / "test_ref.pca_unrelated")


def make_pca_reference(plink2, sites, prefix):
    """PCA basis from a simulated GRCh38 panel, as in src/pca/calculate_pc.sh."""
    rng = random.Random(SEED + 1)
    panel = [f"R{i:03d}" for i in range(1, N_PANEL + 1)]
    with tempfile.TemporaryDirectory() as d:
        vcf = Path(d) / "panel.vcf"
        with open(vcf, "w") as f:
            f.write("##fileformat=VCFv4.2\n"
                    f"##contig=<ID=22,length={B38_LEN['22']}>\n"
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(panel) + "\n")
            for pos, ref, alt, p1, p2 in sorted(sites):
                gts = []
                for i in range(N_PANEL):
                    q = p2 if i >= N_PANEL // 2 else p1
                    gts.append(f"{int(rng.random() < q)}|{int(rng.random() < q)}")
                f.write(f"22\t{pos}\t22:{pos}:{ref}:{alt}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
        run = lambda *a: subprocess.run([plink2, *a], check=True, capture_output=True)
        run("--vcf", vcf, "--make-pgen", "--out", f"{d}/panel")
        run("--pfile", f"{d}/panel", "--maf", "0.05", "--pca", str(N_PC), "allele-wts",
            "--freq", "counts", "--out", f"{d}/pca")
        for ext in (".acount", ".eigenvec.allele", ".eigenval"):
            shutil.copyfile(f"{d}/pca{ext}", f"{prefix}{ext}")


if __name__ == "__main__":
    main()
