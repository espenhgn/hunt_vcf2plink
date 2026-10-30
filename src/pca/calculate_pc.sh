#!/bin/bash

#SBATCH --account=p697_norment_dev
#SBATCH --job-name=refpca
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=7792M
#SBATCH --cpus-per-task=16

set -Eeuo pipefail

root=/ess/p697/cluster/users/alexeas/projects/gnomad
pfile=${root}/pgen/hgdp1kgp_chr_autosomes
# hgdp_1kg_sample_info.unrelateds.pca_outliers_removed.with_project.tsv is downloaded from gnomad,
# "Super population labels TSV" item in "Super population labels TSV" section.
# The file contains QCed, unrelated individuals with information on population.
sample_info=${root}/hgdp_1kg_sample_info.unrelateds.pca_outliers_removed.with_project.tsv
n_pc=40

# NB! Remove extract_snps variable if not testing!!!
extract_snps="/ess/p697/cluster/users/alexeas/projects/husk_imputation/out/husk.genotyped_snps"

plink2="plink2"
threads=16
memory=$((7792 * ${threads}))

echo "Using plink2: $(which ${plink2}), version: $(${plink2} --version)"
echo "Using ${threads} threads and ${memory}Mb RAM."

# Make a list of sample ids.
sample_ids=${sample_info}.ids_only
tail -n+2 ${sample_info} | cut -f1 > ${sample_ids}

# Get common LD independent variants. 
indep_snps=${pfile}.indep_snps
${plink2} --pfile ${pfile} --keep ${sample_ids} --extract ${extract_snps} \
    --indep-pairwise 500kb 0.1 --maf 0.05 --snps-only \
    --out ${indep_snps} --memory ${memory} --threads ${threads}

# Calculate PCs for unrelated samples.
pca_unrelated=${pfile}.pca_unrelated
${plink2} --pfile ${pfile} --keep ${sample_ids} --extract ${indep_snps}.prune.in \
    --pca ${n_pc} 'allele-wts' --freq counts \
    --out ${pca_unrelated} --memory ${memory} --threads ${threads}

# Project PCA for all samples on unrelated PCA basis.
pca_all=${pfile}.pca_all
${plink2} --pfile ${pfile} --read-freq ${pca_unrelated}.acount \
    --score ${pca_unrelated}.eigenvec.allele 2 5 'header-read' 'no-mean-imputation' 'variance-standardize' \
    --score-col-nums 6-$((${n_pc} + 5)) --threads ${threads} --memory ${memory} --out ${pca_all}

