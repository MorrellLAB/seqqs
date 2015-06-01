#!/bin/env bash

#PBS -l mem=4000mb,nodes=1:ppn=4,walltime=8:00:00 
#PBS -m abe 
#PBS -M pmorrell@umn.edu 
#PBS -q lab

set -e
set -u
set -o pipefail

#   The trimming script runs seqqs, scythe, and sickle, and plots_seqqs.R
#   The script is heavily modified from a Vince Buffalo original
#   Most important modification is the addition of plotting of read data before &
#   after 
TRIM_SCRIPT=${HOME}/Apps/Buffalo/seqqs/wrappers/trim_autoplot.sh

#   Extension on forward read names, to be trimmed by basename
READ_NAMING=_1_sequence.txt.gz

#   Project name
PROJECT=SCN

#   Output directory, currently writing full processed directory to scratch
#   Need a symlink at ${HOME} to a scratch directory
OUTDIR=${HOME}/scratch/SCN/rename/sandbox

#   Not currently using this variable
WORKING=${HOME}/scratch/SCN/rename

#   List of samples to be processed
#   Need to hard code the file path for qsub jobs
SAMPLE_INFO=${HOME}/scratch/SCN/rename/SCN_reads.txt
#   Create a bash array of sample names for files
SAMPLE_NAMES=($(cut -f 1 "$SAMPLE_INFO"))
#   Sort this array to make sure that paired end reads end up paired
SAMPLE_NAMES=($(printf '%s\n' "${SAMPLE_NAMES[@]}" | sort))

#   This will list all sample names read into the array
#echo "${SAMPLE_NAMES[@]}"

#   Take the length of array, number of sets of reads, to ensure pairs are \
#   all present
READ_NUMBER=${#SAMPLE_NAMES[@]}

#   This is a pretty simple test, are there an even number of reads?
#   If not, the script should break
PAIRED_READS=$(( $READ_NUMBER % 2 ))
if [ "$PAIRED_READS" -ne 0 ]
   then
echo "An uneven number of paired end reads was found!"
break
fi

#   Assuming we have paired reads, need the count of the penultimate read
#   Will be used below for the iterator
PENULT=$(( $READ_NUMBER -1 ))

#   Change to program directory
#   This is necessary to call the R script (used for plotting) from the trim_autoplot.sh script
cd ${HOME}/Apps/Buffalo/seqqs/wrappers/

#  Bash arrays are zero indexed, so the sequence has to be 0, 2, 4 for forward reads
#  For reverse reads it will be 1 3 5 to the last element
for i in $(seq 0 2 $PENULT)
    do
SAMPLE_NAME=$(basename ${SAMPLE_NAMES[i]} ${READ_NAMING})

        #   Path to forward reads
    FORWARD_READS=${WORKING}/${SAMPLE_NAMES[i]}
        #   Below should iterate over reverse read positions
        #   Don't have an explicit test that reads are paired!
for j in $(seq 1 2 $READ_NUMBER)
    do
        #   Path to reverse reads
    REVERSE_READS=${WORKING}/${SAMPLE_NAMES[j]}

# Run the script and write output to a directory defined by project and sample name
${TRIM_SCRIPT} ${SAMPLE_NAME} ${FORWARD_READS} ${REVERSE_READS} ${OUTDIR}/${PROJECT}/${SAMPLE_NAME}/
done
done
