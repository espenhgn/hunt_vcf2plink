#!/bin/bash
# Convert genotypes from VCF format to PLINK2 binary (pgen) file format
# run "mamba activate biobanknorge" before invoking

export THREADS=8
export MEMPERTHREAD=1948
export MEMORY=$(( ${MEMPERTHREAD} * ${THREADS} ))

set -Eeuo pipefail

export APPTAINER_TMPDIR=$TMPDIR

export ROOTDIR="/home/ubuntu/shared"
export SOURCEDIR="${ROOTDIR}/data/vcf"
export DESTDIR="${ROOTDIR}/data/plink2.GRCh37"
export SIF="${ROOTDIR}/github/comorment/containers/containers/1.16.0/python3.sif"
export PLINK2="apptainer exec --bind=/home/ubuntu:/home/ubuntu ${SIF} plink2"

mkdir -p $DESTDIR

# Function to process a single chromosome
process_chromosome() {
    local CHR=$1
    export PREFIX="CHR${CHR}_PID107822"
    export VCF=${SOURCEDIR}/${PREFIX}.vcf.gz
    export PFILE=${DESTDIR}/${PREFIX}

    if [[ "$CHR" == "X" ]]; then
        ${PLINK2} \
            --vcf ${VCF} dosage=DS \
            --make-pgen \
            --maf 0.005 \
            --exclude-if-info 'INFO<0.3' \
            --geno 0.05 \
            --split-par b37 \
            --hwe 1E-50 0.001 midp \
            --update-sex \
            --out ${PFILE} \
            --threads ${THREADS} \
            --memory ${MEMORY} | tee /dev/stderr
    else
        ${PLINK2} \
            --vcf ${VCF} dosage=DS \
            --make-pgen \
            --maf 0.005 \
            --exclude-if-info 'INFO<0.3' \
            --geno 0.05 \
            --split-par b37 \
            --hwe 1E-50 0.001 midp \
            --out ${PFILE} \
            --threads ${THREADS} \
            --memory ${MEMORY} | tee /dev/stderr
    fi

}

# Export the function and necessary variables for GNU parallel
export -f process_chromosome

# Use GNU parallel to process chromosomes concurrently
# seq 3 22 | parallel -j 4 process_chromosome

# Use GNU parallel to process specific chromosomes concurrently
parallel -j 4 process_chromosome ::: 1 7 12