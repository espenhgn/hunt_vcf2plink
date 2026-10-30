#!/bin/bash
# Align SNP alleles to the GRCh38 reference sequence and rename variants to
# "chr:pos:ref:alt", as needed e.g. for PCA scoring. Adapted from src/align_ids_toGRCh38.sh.
# Indels, non-ACGT, strand ambiguous and monomorphic variants are dropped.
# On chrX the PAR is merged into X for alignment (no PAR contig in the FASTA)
# and split off again (b38) in the final file set.
#
# usage: align.sh PLINK2 REF_FASTA IN_PREFIX OUT_PREFIX CHR THREADS MEMORY WORKDIR SUMMARY

set -Eeuo pipefail

PLINK2=$1
REF_FASTA=$2
PFILE=$3
OUT=$4
CHR=$5
THREADS=$6
MEMORY=$7
WORK=$8
SUMMARY=$9

mkdir -p "${WORK}"

merge_par=""
split_par=""
if [[ "${CHR}" == "X" ]]; then
    merge_par="--merge-par"
    split_par="--split-par b38"
fi

# Perform basic QC: remove indels, non-atgc, strand ambiguous and monomorphic variants and drop non-standard chromosomes.
strand_ambig=${WORK}/strand_ambig.snps
awk 'BEGIN{cpl["A"]="T";cpl["T"]="A";cpl["G"]="C";cpl["C"]="G"} !/^#/ {if(cpl[$4]==$5) print($3)}' "${PFILE}.pvar" > "${strand_ambig}"
pfile_qc=${WORK}/qc
${PLINK2} --pfile "${PFILE}" ${merge_par} \
    --make-pgen 'pvar-cols=' \
    --out "${pfile_qc}" \
    --exclude "${strand_ambig}" \
    --chr '1-26' \
    --snps-only 'just-acgt' \
    --rm-dup 'exclude-mismatch' \
    --mac 1 \
    --threads "${THREADS}" \
    --memory "${MEMORY}"

# Align variant alleles to the reference genome sequence.
# This is assumed to be done after the qc step, i.e. dataset contains no indels, no non-atgc, no strand ambiguous and no monomorphic variants.
# --ref-from-fa orders ref and alt alleles so that ref allele matches allele in the reference sequence.
# If neither of two alleles fits the reference, it adds PR (provisional reference) tag to the INFO column.
# N.b. with --ref-from-fa only reordering but no strand flip is attempted, therefore strand flip is implemented manually.
pfile_align=${WORK}/align
${PLINK2} --pfile "${pfile_qc}" \
    --make-pgen 'pvar-cols=+info' \
    --out "${pfile_align}" \
    --ref-from-fa 'force' \
    --fa "${REF_FASTA}" \
    --threads "${THREADS}" \
    --memory "${MEMORY}"

# Flip the strand of variants which contain PR tag in the INFO column, i.e. those which were
# not aligned to the reference sequence in the previous command.
pfile_align_flip=${pfile_align}.flip
awk 'BEGIN{OFS="\t"; cpl["A"]="T";cpl["T"]="A";cpl["G"]="C";cpl["C"]="G"} {if($0!~/^#/ && $6=="PR") {$4=cpl[$4]; $5=cpl[$5]} print($0)}' "${pfile_align}.pvar" > "${pfile_align_flip}.pvar"

# Align variant alleles (which are now flipped) to reference sequence again.
pfile_align_flip_align=${pfile_align_flip}.align
${PLINK2} --pvar "${pfile_align_flip}.pvar" \
    --psam "${pfile_align}.psam" \
    --pgen "${pfile_align}.pgen" \
    --make-pgen 'pvar-cols=+info' \
    --out "${pfile_align_flip_align}" \
    --ref-from-fa 'force' \
    --fa "${REF_FASTA}" \
    --threads "${THREADS}" \
    --memory "${MEMORY}"

# Make input file for --update-name taking only variants without PR tag in the INFO column,
# i.e. those with alleles correctly aligned to the reference sequence.
# The input file for --update-name contains two columns:
#   (1) original variant ids;
#   (2) new "chr:bp:ref:alt", where ref allele is the allele in the reference sequence.
update_name_snps=${WORK}/update_name.snps
awk 'BEGIN{OFS="\t"} {if($0!~/^#/ && $6!="PR") print($3,$1":"$2":"$4":"$5)}' "${pfile_align_flip_align}.pvar" > "${update_name_snps}"
# Make a list of SNPs to exclude, i.e. those for which alleles were not aligned to the
# reference sequence neither on the original nor on the flipped strand.
exclude_snps=${WORK}/exclude.snps
awk '{if($0!~/^#/ && $6=="PR") print($3)}' "${pfile_align_flip_align}.pvar" > "${exclude_snps}"

# Make pgen file with "chr:bp:ref:alt" variant IDs keeping only variants which were aligned
# to the reference sequence.
# This step relies on the order of operations in plink2:
# https://www.cog-genomics.org/plink/2.0/order.
# First goes --update-name, then --exclude and finally --rm-dup.
${PLINK2} --pfile "${pfile_align_flip_align}" \
    --make-pgen \
    --out "${OUT}" \
    --update-name "${update_name_snps}" \
    --exclude "${exclude_snps}" \
    --rm-dup 'force-first' 'list' \
    ${split_par} \
    --threads "${THREADS}" \
    --memory "${MEMORY}"

{
    printf 'chr\t%s\n' "${CHR}"
    printf 'align_input\t%s\n' "$(grep -vc '^#' "${PFILE}.pvar")"
    printf 'align_after_qc\t%s\n' "$(grep -vc '^#' "${pfile_qc}.pvar")"
    printf 'align_flipped\t%s\n' "$(awk '!/^#/ && $6=="PR" {n++} END {print n+0}' "${pfile_align}.pvar")"
    printf 'align_dropped_unaligned\t%s\n' "$(awk 'END {print NR}' "${exclude_snps}")"
    printf 'align_output\t%s\n' "$(grep -vc '^#' "${OUT}.pvar")"
} > "${SUMMARY}"

rm -rf "${WORK:?}"
