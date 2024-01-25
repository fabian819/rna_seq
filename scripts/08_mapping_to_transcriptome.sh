#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --partition=pall
#SBATCH --mem-per-cpu=8000M

FDIR="./data/fastq/"
ADIR="/data/users/fgribi/rna_seq/data/annotations/GRCh38_APPRIS_CDS_18/GRCh38_APPRIS_CDS_18"
BDIR="/data/users/fgribi/rna_seq/data/bam/"

cd $FDIR

module load UHTS/Aligner/bowtie/1.2.0

for x in $(ls -d *RNA.fastq); 
do echo ${x};
bowtie \
 -t \
 -p 4 \
 -v 1 \
 -m 1 \
 --best \
 --strata \
 --norc \
 $ADIR \
 -q ${x} \
 -S ${BDIR}/$(basename ${x} .fastq)_GRCh38_APPRIS_CDS.sam 2>$(basename ${x} .fastq)_GRCh38_APPRIS_CDS_log.txt; done
