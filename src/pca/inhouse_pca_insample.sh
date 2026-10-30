#!/bin/bash

# compute PC scores using eigenvectors/values computed in independent sample for HUNT sample
# run "mamba activate biobanknorge" before invoking

set -o errexit

export APPTAINER_TMPDIR=$TMPDIR

export ROOTDIR="/home/ubuntu/shared"                                                                                                                      
export INPUTDIR="${ROOTDIR}/data/GRCh38.aligned.merged"
export OUTPUTDIR="${INPUTDIR}.pca"
export PPREFIX="CHR1_22_PID107822"
export SIF="${ROOTDIR}/github/comorment/containers/containers/1.11.0/python3.sif"                                                               
export PLINK2="apptainer exec --bind=/home/:/home/ ${SIF} plink2"

threads=16
mem_per_cpu=1948

echo "Using $(${PLINK2} --version)"
memory=$((threads * mem_per_cpu))
echo "Using ${threads} threads and ${memory}Mb RAM."

# Input arguments.
pfile=${INPUTDIR}/${PPREFIX}
# bfile=${pfile}
# The following files are created using "calculate_pc.sh" script for reference processing.
# cf. w. Alexey on how this was originally done on p697 for Husk sample
ref_root="${ROOTDIR}/data/scripts/pca/pgen"
ref_pfile=${ref_root}/hgdp1kgp_chr_autosomes
ref_acount=${ref_pfile}.pca_unrelated.acount
ref_eigenvec=${ref_pfile}.pca_unrelated.eigenvec.allele
n_pc=40 # should be less or equal to the number of PCs calculated for the reference panel

# Project PCA for all samples on unrelated PCA basis
mkdir -p $OUTPUTDIR
$PLINK2 --pfile ${pfile} --read-freq ${ref_acount} \
    --score ${ref_eigenvec} 2 5 'header-read' 'no-mean-imputation' 'variance-standardize' \
    --score-col-nums 6-$((${n_pc} + 5)) \
    --threads ${threads} \
    --memory ${memory} \
    --out ${OUTPUTDIR}/${PPREFIX}.pca_all
    