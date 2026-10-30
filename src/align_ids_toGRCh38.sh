#!/bin/bash
# Align genotypes to known reference panel, needed e.g., by PCA scoring scripts
# run "mamba activate biobanknorge" before invoking

set -Eeuo pipefail

export THREADS=8
export MEMPERTHREAD=1948
export MEMORY=$(( ${MEMPERTHREAD} * ${THREADS} ))

export APPTAINER_TMPDIR=$TMPDIR

export ROOTDIR="/home/ubuntu/shared"
export SOURCEDIR="${ROOTDIR}/data/plink2.GRCh38"
export DESTDIR="${SOURCEDIR}.aligned"
export SIF="${ROOTDIR}/github/comorment/containers/containers/1.11.0/python3.sif"
export PLINK2="apptainer exec --bind=/home/:/home/ ${SIF} plink2"

mkdir -p $DESTDIR

# export REF_FASTA="${ROOTDIR}/reference/data/hs37d5.fa.zst"  # build 37
export REF_FASTA="${ROOTDIR}/data/reference/data/GRCh38_full_analysis_set_plus_decoy_hla.fa.zst" # build 38

# Function to process a single chromosome
process_chromosome() {
    local CHR=$1
    export PREFIX="CHR${CHR}_PID107822"
    export PFILE=${SOURCEDIR}/${PREFIX}

    ${PLINK2} --version

    # Perform basic QC: remove indels, non-atgc, strand ambiguous and monomorphic variants and drop non-standard chromosomes.
    strand_ambig=${TMPDIR}/tmp_chr${CHR}.strand_ambig.snps
    awk 'BEGIN{cpl["A"]="T";cpl["T"]="A";cpl["G"]="C";cpl["C"]="G"} {if(cpl[$4]==$5) print($3)}' ${PFILE}.pvar > ${strand_ambig}
    pfile_qc=${TMPDIR}/tmp_chr${CHR}.qc
    ${PLINK2} --pfile ${PFILE} \
        --make-pgen 'pvar-cols=' 'psam-cols=' \
        --out ${pfile_qc} \
        --exclude ${strand_ambig} \
        --chr '1-26' \
        --snps-only 'just-acgt' \
        --rm-dup 'exclude-mismatch' \
        --mac 1 \
        --threads ${THREADS} \
        --memory ${MEMORY}

    # Align variant alleles to the reference genome sequence.
    # This is assumed to be done after the qc step, i.e. dataset contains no indels, no non-atgc, no strand ambiguous and no monomorphic variants. 
    # --ref-from-fa orders ref and alt alleles so that ref allele matches allele in the reference sequence.
    # If neither of two alleles fits the reference, it adds PR (provisional reference) tag to the INFO column.
    # N.b. with --ref-from-fa only reordering but no strand flip is attempted, therefore strnd flip is implemented manually.
    pfile_align=${TMPDIR}/tmp_chr${CHR}.align
    ${PLINK2} --pfile ${pfile_qc} \
        --make-pgen 'pvar-cols=+info' \
        --out ${pfile_align} \
        --ref-from-fa 'force' \
        --fa ${REF_FASTA} \
        --threads ${THREADS} \
        --memory ${MEMORY}

    # Flip the strand of variants which contain PR tag in the INFO column, i.e. those which were
    # not aligned to the reference sequence in the previous command.
    pfile_align_flip=${pfile_align}.flip
    awk 'BEGIN{OFS="\t"; cpl["A"]="T";cpl["T"]="A";cpl["G"]="C";cpl["C"]="G"} {if($0!~/^#/ && $6=="PR") {$4=cpl[$4]; $5=cpl[$5]} print($0)}' ${pfile_align}.pvar > ${pfile_align_flip}.pvar

    # Align variant alleles (which are now flipped) to reference sequence again.
    pfile_align_flip_align=${pfile_align_flip}.align
    ${PLINK2} --pvar ${pfile_align_flip}.pvar \
        --psam ${pfile_align}.psam \
        --pgen ${pfile_align}.pgen \
        --make-pgen 'pvar-cols=+info' \
        --out ${pfile_align_flip_align} \
        --ref-from-fa 'force' \
        --fa ${REF_FASTA} \
        --threads ${THREADS} \
        --memory ${MEMORY} 

    # Make input file for --update-name taking only variants without PR tag in the INFO column,
    # i.e. those with alleles correctly aligned to the reference sequence.
    # The input file for --update-name contains two columns:
    #   (1) original variant ids;
    #   (2) new "chr:bp:ref:alt", where ref allele is the allele in the reference sequence.
    update_name_snps=${TMPDIR}/tmp_chr${CHR}.update_name.snps
    awk 'BEGIN{OFS="\t"} {if($0!~/^#/ && $6!="PR") print($3,$1":"$2":"$4":"$5)}' ${pfile_align_flip_align}.pvar > ${update_name_snps}
    # Make a list of SNPs to exclude, i.e. those for which alleles were not aligned to the
    # reference sequence neither on the original nor on the flipped strand. 
    exclude_snps=${TMPDIR}/tmp_chr${CHR}.exclude.snps
    touch ${exclude_snps} # ensure the file exists, even if there are no SNPs with allele mismatch.
    awk '{if($0!~/^#/ && $6=="PR") print($3)}' ${pfile_align_flip_align}.pvar > ${exclude_snps}

    # Make pgen file with "chr:bp:ref:alt" variant IDs keeping only variants which were aligned
    # to the reference sequence.
    # This step relies on the order of operations in plink2:
    # https://www.cog-genomics.org/plink/2.0/order.
    # First goes --update-name, then --exclude and finally --rm-dup.
    aligned_pfile=${DESTDIR}/${PREFIX}
    ${PLINK2} --pfile ${pfile_align_flip_align} \
        --make-pgen \
        --out ${aligned_pfile} \
        --update-name ${update_name_snps} \
        --exclude ${exclude_snps} \
        --rm-dup 'force-first' 'list' \
        --threads ${THREADS} \
        --memory ${MEMORY}

    # clean up tmp files
    rm -f ${strand_ambig} ${pfile_qc}.* ${pfile_align}.* ${pfile_align_flip}.* ${update_name_snps} ${exclude_snps} || true

}

# Export the function and necessary variables for GNU parallel
export -f process_chromosome

# Use GNU parallel to process chromosomes concurrently
seq 1 22 | parallel -j 4 process_chromosome
# parallel -j 4 process_chromosome ::: 19
