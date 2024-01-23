#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --partition=pall

module add UHTS/Analysis/sratoolkit/2.10.7

TDIR="./data"
NAMES="SRR9596295 SRR9596296 SRR9596300 SRR9596310 SRR9596303 SRR9596304"

if [ ! -d $TDIR ]; then
    mkdir $TDIR
fi

cd $TDIR

#Fetch and convert the desired data
for i in $NAMES;
do
    prefetch ${i}; 
    fastq-dump --gzip ${i};
done

rm */*.sra
rmdir SRR*