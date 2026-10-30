#!/bin/sh
# merge per-chromosome files to both PLINK1 and PLINK2 format
# run "mamba activate biobanknorge" before invoking

# module purge
set -o errexit

export APPTAINER_TMPDIR=$TMPDIR

export ROOTDIR="/home/ubuntu/shared"
# export BUILD="GRCh37"                                                
export BUILD="GRCh38"                                                                           
export SOURCEDIR="${ROOTDIR}/data"
export INPUTDIR="${SOURCEDIR}/plink2.${BUILD}.aligned"
export OUTPUTDIR="${SOURCEDIR}/${BUILD}.aligned.merged"
export SIF="${ROOTDIR}/github/comorment/containers/containers/1.11.0/python3.sif"
export PLINK2="apptainer exec --bind=/home/:/home/ ${SIF} plink2"

export MERGED='CHR1_22_PID107822'

mkdir -p $OUTPUTDIR
rm -f ${OUTPUTDIR}/mergelist.txt

for i in $(seq 1 22)
do
    export PREFIX="CHR${i}_PID107822" 
    echo "${INPUTDIR}/$PREFIX" >> ${OUTPUTDIR}/mergelist.txt
done

$PLINK2 --pmerge-list ${OUTPUTDIR}/mergelist.txt pfile \
    --threads 16 \
    --memory 32000 \
    --make-bed \
    --out ${OUTPUTDIR}/${MERGED}

