#!/bin/sh
# lift per-chromosome PLINK2 filesets to build GRCh38
# run "mamba activate biobanknorge" before invoking

set -o errexit

export APPTAINER_TMPDIR=$TMPDIR

export ROOTDIR="/home/ubuntu/shared"                                                                                                              
export SOURCEDIR="${ROOTDIR}/data"                                                                                     
export INPUTDIR="${SOURCEDIR}/plink2.GRCh37"
export OUTPUTDIR="${SOURCEDIR}/plink2.GRCh38"
export SIF="${ROOTDIR}/github/comorment/containers/containers/1.11.0/python3.sif"
# export PLINK="apptainer exec --bind=/home:/home ${SIF} plink"
# export PLINK2="apptainer exec --bind=/home:/home ${SIF} plink2"
export PLINK2="/home/conda/envs/biobanknorge/bin/plink2"
export THREADS=8
export MEMPERTHREAD=1948
export MEMORY=$(( ${MEMPERTHREAD} * ${THREADS} ))

mkdir -p $OUTPUTDIR

process_chromosome() {
    local CHR=$1
    # input 
    export PREFIX="CHR${CHR}_PID107822"
    export PFILE="${INPUTDIR}/${PREFIX}"

    # convert to plink1 format assumed by tripleliftover
    $PLINK2 --pfile $PFILE  \
        --threads $THREADS \
        --memory $MEMORY \
        --rm-dup 'exclude-mismatch' \
        --make-just-bim \
        --out $TMPDIR/$PREFIX 

    # Run triple-liftOver for your bim file with hg19 as base build and hg38 as target build:
    perl ${ROOTDIR}/github/GraceSheng/triple-liftOver/tripleliftover_v133.pl \
        --bim $TMPDIR/$PREFIX.bim \
        --base hg19 \
        --target hg38 \
        --outprefix ${TMPDIR}/${PREFIX}

    # use the variant name in the output to exclude any variant that is not in the 
    # triple-liftOver output (these variants either failed the lift over or were lifted to a 
    # different chromosome). It can be executed with a command in Linux followed by a PLINK command:
    awk '{if (NR>1) print $1}' "${TMPDIR}/${PREFIX}.invr.txt" > $TMPDIR/${PREFIX}.variant.name.in.liftOver.output

    # Use the variant name and pos_[Hg38] from the output to update the positions in the hg19 PLINK files. 
    # This step generates the hg38 PLINK file with genome forward strand alleles 
    # (with the inverted site issue between hg19 and hg38):
    awk '{if (NR>1) print $1 "\t" $3}' "${TMPDIR}/${PREFIX}.invr.txt" > $TMPDIR/${PREFIX}.variant.name.hg38.position

    # Parse out the “inverted” variant names from triple-liftOver output, limit them to SNVs only 
    # (as it may include indels, exclusion of the indels is not shown in the example). 
    # Execute the strand flip command in PLINK for these inverted SNVs. 
    # This will generate the final hg38 PLINK dataset for imputation (fixes the inverted site issue 
    # between hg19 and hg38):
    awk '{if (NR>1 && $5 == "inverted") print $1}' "${TMPDIR}/${PREFIX}.invr.txt" > $TMPDIR/${PREFIX}.inverted.SNV.name
    
    # Combine all flags in a single pass, not sequential as suggested on github.com/GraceSheng/triple-liftOver 
    # (md5 checksums/results are indistinguashible)
    $PLINK2 --pfile "${PFILE}" \
        --extract $TMPDIR/${PREFIX}.variant.name.in.liftOver.output \
        --update-map $TMPDIR/${PREFIX}.variant.name.hg38.position \
        --flip $TMPDIR/${PREFIX}.inverted.SNV.name \
        --sort-vars \
        --threads $THREADS \
        --memory $MEMORY \
        --make-pgen \
        --out ${OUTPUTDIR}/${PREFIX}

    # clean up
    rm $TMPDIR/${PREFIX}.*

}

# Export the function and necessary variables for GNU parallel
export -f process_chromosome

# Use GNU parallel to process chromosomes concurrently
seq 1 22 | parallel -j 4 process_chromosome
# parallel -j 4 process_chromosome ::: 19
# process_chromosome 19