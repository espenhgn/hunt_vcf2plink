#!/bin/sh
#SBATCH --account=p2938_tsd
#SBATCH --error=/ess/p2938/data/durable/tromsoestudy/genetics/dementia/logs/merge_chr.txt
#SBATCH --output=/ess/p2938/data/durable/tromsoestudy/genetics/dementia/logs/merge_chr.txt
#SBATCH --open-mode=truncate
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4000M
# merge per-chromosome files to both PLINK1 and PLINK2 format

module purge
set -o errexit

export ROOTDIR="/ess/p2938/data/durable"                                                                                                                      
export SOURCEDIR="${ROOTDIR}/tromsoestudy/genetics/dementia/TU_imputed_8030_00676"
export BUILD="GRCh37"                                                                              
export INPUTDIR="${SOURCEDIR}.plink2.${BUILD}.aligned"
export OUTPUTDIR="${SOURCEDIR}.${BUILD}.aligned.merged"
export SIF="${ROOTDIR}/users/espehage/github/comorment/containers/containers/1.9.1/python3.sif"                                                               
export PLINK2="apptainer exec --bind=/ess:/ess ${SIF} plink2"

export MERGED=gwas_data_merged

mkdir -p $OUTPUTDIR
rm -f ${OUTPUTDIR}/mergelist.txt

for i in $(seq 1 22)
do
    export PREFIX="gwas_data_chr${i}.renamed"
    echo "${INPUTDIR}/$PREFIX" >> ${OUTPUTDIR}/mergelist.txt
done

$PLINK2 --pmerge-list ${OUTPUTDIR}/mergelist.txt pfile \
    --threads 8 \
    --memory 32000 \
    --make-bed \
    --out ${OUTPUTDIR}/${MERGED}

