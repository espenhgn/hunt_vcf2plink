# hunt_vcf2plink

Snakemake workflow that converts the per-chromosome HUNT VCF files (GRCh37, imputed
dosages) to PLINK file sets in GRCh38, and projects samples onto a reference PCA basis
(HGDP + 1000 Genomes, PC1–40). It runs in a container in HUNT Cloud.

## Steps

Per chromosome (1–22, X):

1. **VCF → pgen (GRCh37)**: dosages from `DS`; QC: `INFO >= 0.3`, `MAF >= 0.005`,
   `--geno 0.05`, `--hwe 1E-50 0.001 midp`. Sex comes from the covariates file
   (`Sex`: 1 = male, 2 = female). The PAR is split off chrX (b37) so that male PAR
   genotypes are diploid in QC.
2. **Liftover to GRCh38** with [triple-liftOver](https://github.com/GraceSheng/triple-liftOver).
   Variants that don't lift, or lift to another chromosome, are dropped. SNVs in
   regions inverted between builds are strand flipped. Inverted indels are dropped,
   as are variants with duplicated IDs. On chrX the PAR is merged into X for the
   liftover and split off again afterwards (b38).
3. **GRCh38, all variants**: IDs are set to `chr:pos:ref:alt` in GRCh38 coordinates.
   Alleles are as in the input, not aligned to the reference.
4. **GRCh38, aligned**: SNPs only, with REF/ALT aligned (and strand flipped where
   needed) to the GRCh38 reference sequence, and IDs `chr:pos:ref:alt`. Indels,
   non-ACGT, strand-ambiguous and monomorphic variants are dropped.

Then:

- **Merge** the aligned per-chromosome file sets.
- **PCA projection** onto the HGDP+1KG basis in [src/pca/pgen](src/pca/pgen)
   (autosomal SNPs, GRCh38). The workflow fails if fewer than `PCA_MIN_OVERLAP` of
   the reference variants are found in the data.
- **Report** with variant counts per step and chromosome.

The workflow is adapted from the original scripts in [src/](src/). `.pgen` file sets
keep dosages. `.bed` file sets have hard calls (`QC.hard_call_threshold`) and
numeric chromosome codes (23 = X, 25 = PAR) for PLINK 1 tools. In `.pvar` files the
PAR is coded `PAR1`/`PAR2`, as in PLINK 2.0 a7. PAR variant IDs start with `X:`.

## Outputs (`OUT_DIR`)

```text
GRCh38/CHR{chr}_PID107822.{pgen,pvar,psam,bed,bim,fam}         all lifted variants, per chromosome
GRCh38.aligned/CHR{chr}_PID107822.{pgen,pvar,psam}             aligned SNPs, per chromosome
GRCh38.aligned/CHR1_22X_PID107822.{pgen,pvar,psam,bed,bim,fam} aligned SNPs, merged
pca/CHR1_22X_PID107822.pcs.tsv                                 IID, PC1..PC40
pca/CHR1_22X_PID107822.sscore                                  plink2 --score output
qc/report.tsv                                                  counts per step
logs/                                                          per-rule logs
```

## Running in HUNT Cloud

1. Get the container image. Either pull it (if the lab allows outbound access to ghcr.io):

   ```bash
   apptainer pull hunt_vcf2plink_0.1.0.sif docker://ghcr.io/espenhgn/hunt_vcf2plink:0.1.0
   ```

   or build the `.sif` outside the lab and import it through the HUNT Cloud file
   import.
2. Edit [config.yaml](config.yaml) (paths, `CHROMS`, resources) and
   [profiles/hunt/config.yaml](profiles/hunt/config.yaml) (cores, memory of the VM).
3. Run from the repository checkout:

   ```bash
   cd /home/ubuntu/shared/github/espenhgn/hunt_vcf2plink
   apptainer exec --bind /home/ubuntu/shared hunt_vcf2plink_0.1.0.sif \
       snakemake --profile profiles/hunt -n   # dry run first
   apptainer exec --bind /home/ubuntu/shared hunt_vcf2plink_0.1.0.sif \
       snakemake --profile profiles/hunt
   ```

   To try a single chromosome first: `--config 'CHROMS=[22]'`.

Snakemake can also run on the host and start the container per job, using
`--software-deployment-method apptainer` (the image is set by `DOCKER` in the config).

## Requirements

All tools are in the container ([docker/Dockerfile](docker/Dockerfile),
[environment.yml](environment.yml)): Snakemake 9.27, PLINK 2.0 a7.8,
triple-liftOver (commit `f781ed6`) with UCSC liftOver 482 and the
`hg19ToHg38` chain file. The GRCh38 reference FASTA
(`GRCh38_full_analysis_set_plus_decoy_hla.fa.zst`) is not in the image; set
`REF_FASTA_GRCH38`.

UCSC liftOver and the chain files are free for non-commercial use only; see the
[UCSC license](https://genome.ucsc.edu/license/).

## Testing

[test/](test/) holds a small synthetic data set (100 samples, chr22 and chrX; no real
data) made by [test/make_test_data.py](test/make_test_data.py), and a config for it:

```bash
docker build --platform linux/amd64 -f docker/Dockerfile -t hunt_vcf2plink:test .
docker run --rm -v "$PWD:/work" -w /work --user "$(id -u):$(id -g)" -e HOME=/tmp \
    hunt_vcf2plink:test snakemake -s workflow/Snakefile --configfile test/config.yaml --cores 2
```

CI ([.github/workflows/docker.yml](.github/workflows/docker.yml)) builds the image,
runs this test, and pushes the image to ghcr.io on pushes to `main` and on `v*` tags.

## Input VCF files

```bash
$ ls -l ~/shared/data/vcf/
total 540454964
-rw-rw-r-- 1 ubuntu ubuntu 25626482793 Jul  8  2025 CHR10_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 25191496744 Jul  8  2025 CHR11_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 24307527149 Jul 10  2025 CHR12_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 18415273991 Jul  8  2025 CHR13_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 16999844704 Jul  8  2025 CHR14_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 15406032393 Jul  8  2025 CHR15_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 17346277497 Jul  8  2025 CHR16_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 15321765020 Jul  8  2025 CHR17_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 15158074313 Jul  8  2025 CHR18_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 13058352845 Jul  8  2025 CHR19_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 39685604883 Jul  8  2025 CHR1_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 12151436979 Jul  8  2025 CHR20_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu  7530562842 Jul  8  2025 CHR21_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu  8055908367 Jul  8  2025 CHR22_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 42631353016 Jul  8  2025 CHR2_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 36177923160 Jul  8  2025 CHR3_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 36579476123 Jul  8  2025 CHR4_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 32556604281 Jul  7  2025 CHR5_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 35447744795 Jul  8  2025 CHR6_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 30186180253 Jul  8  2025 CHR7_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 28171206523 Jul  8  2025 CHR8_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 22335290871 Jul  8  2025 CHR9_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 35085262552 Jul  9  2025 CHRX_PID107822.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu        1301 Jul  7  2025 md5_vcf_PID107822.txt
```

## Citations

- PLINK 2.0: Chang et al. (2015) GigaScience 4:7.
- triple-liftOver: Sheng et al. (2022) bioRxiv, doi:10.1101/2022.02.19.481174.
- HGDP + 1KG reference: Koenig et al. (2024) Genome Research (gnomAD HGDP+1KG callset).
