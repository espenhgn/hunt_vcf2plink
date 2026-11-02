#!/bin/bash
# Lift a GRCh37 pgen file set to GRCh38 with UCSC liftOver, using the triple-liftOver
# method to find sites in inverted regions (liftover_bim.py). Adapted from
# src/plink2_liftover.sh.
#
# Variants that fail to lift, or lift to another chromosome, are dropped. SNVs in
# regions inverted between builds are strand flipped; inverted non-SNVs are dropped,
# as flipping does not fix them. On chrX the PAR is merged back into X first, so that
# PAR variants are lifted and position-sorted together with X.
#
# usage: liftover.sh PLINK2 LIFTOVER CHAIN IN_PREFIX OUT_PREFIX CHR THREADS MEMORY WORKDIR SUMMARY

set -Eeuo pipefail

PLINK2=$1
LIFTOVER=$2
CHAIN=$3
IN=$4
OUT=$5
CHR=$6
THREADS=$7
MEMORY=$8
WORK=$9
SUMMARY=${10}
SCRIPTS=$(dirname "$0")

mkdir -p "${WORK}"

merge_par=""
if [[ "${CHR}" == "X" ]]; then
    merge_par="--merge-par"
fi

# Drop all variants with duplicated IDs: the liftover output is keyed by ID, and
# --update-map fails on duplicates even when --rm-dup is given in the same run.
awk '!/^#/ {n[$3]++} END {for (id in n) if (n[id] > 1) print id}' "${IN}.pvar" > "${WORK}/duplicated_ids.txt"
n_dup_ids=$(awk 'END {print NR}' "${WORK}/duplicated_ids.txt")

${PLINK2} --pfile "${IN}" ${merge_par} \
    --exclude "${WORK}/duplicated_ids.txt" \
    --make-just-bim \
    --threads "${THREADS}" \
    --memory "${MEMORY}" \
    --out "${WORK}/in"

# columns: snpid chr pos_target pos_base category (lifted/inverted), with header
invr="${WORK}/lift.invr.txt"
python "${SCRIPTS}/liftover_bim.py" \
    --bim "${WORK}/in.bim" \
    --liftover "${LIFTOVER}" \
    --chain "${CHAIN}" \
    --work-dir "${WORK}/liftover" \
    --out "${invr}"

awk -v w="${WORK}" '
    NR == FNR {snv[$2] = ($5 ~ /^[ACGT]$/ && $6 ~ /^[ACGT]$/); next}
    FNR == 1 {next}
    $5 == "inverted" && !snv[$1] {print $1 > (w "/inverted_nonsnv.txt"); next}
    {print $1 > (w "/extract.txt"); print $1 "\t" $3 > (w "/update_map.txt")}
    $5 == "inverted" {print $1 > (w "/flip.txt")}' \
    "${WORK}/in.bim" "${invr}"
touch "${WORK}/extract.txt" "${WORK}/update_map.txt" "${WORK}/flip.txt" "${WORK}/inverted_nonsnv.txt"

n_in=$(awk 'END {print NR}' "${WORK}/in.bim")
n_lifted=$(awk 'END {print NR}' "${WORK}/extract.txt")
if [[ "${n_lifted}" -eq 0 ]]; then
    echo "ERROR: no variants lifted over to hg38 for chr${CHR}" >&2
    exit 1
fi

# combine all flags in a single pass (results identical to the sequential steps
# suggested for triple-liftOver)
${PLINK2} --pfile "${IN}" ${merge_par} \
    --extract "${WORK}/extract.txt" \
    --update-map "${WORK}/update_map.txt" \
    --flip "${WORK}/flip.txt" \
    --sort-vars \
    --make-pgen \
    --threads "${THREADS}" \
    --memory "${MEMORY}" \
    --out "${OUT}"

{
    printf 'chr\t%s\n' "${CHR}"
    printf 'liftover_dropped_duplicate_ids\t%s\n' "${n_dup_ids}"
    printf 'liftover_input\t%s\n' "${n_in}"
    printf 'liftover_lifted\t%s\n' "${n_lifted}"
    printf 'liftover_flipped\t%s\n' "$(awk 'END {print NR}' "${WORK}/flip.txt")"
    printf 'liftover_dropped_inverted_nonsnv\t%s\n' "$(awk 'END {print NR}' "${WORK}/inverted_nonsnv.txt")"
    printf 'liftover_dropped_unmapped\t%s\n' "$(( n_in - n_lifted - $(awk 'END {print NR}' "${WORK}/inverted_nonsnv.txt") ))"
} > "${SUMMARY}"

rm -rf "${WORK:?}"
