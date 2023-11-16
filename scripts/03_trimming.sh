#!/usr/bin/env bash
#SBATCH --cpus-per-task=4 
#SBATCH --time=03:00:00
#SBATCH --partition=pall

module load UHTS/Quality_control/cutadapt/2.5

FDIR="./data/fastq"
FILES="RPF_KO_Rep1 RPF_KO_Rep2 RPF_WT_Rep1 RPF_WT_Rep2"

cd $FDIR

for i in $FILES; do
    #Clip the 3' adapter sequence
    cutadapt \
    -j 4 \
    -a NNNNCTGTAGGCACCATCAAT \
    -q 25 \
    --minimum-length 25 \
    --discard-untrimmed \
    -o "${i}_clpd.fastq.gz" \
    "${i}.fastq.gz"
done