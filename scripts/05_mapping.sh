#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --partition=pall
#SBATCH --mem-per-cpu=8000M

FDIR="./data/fastq/"
ADIR="/data/users/fgribi/rna_seq/data/annotations/GRCh38.dna.primary_assembly/GRCh38.dna.primary_assembly"
BDIR="/data/users/fgribi/rna_seq/data/bam/"

cd $FDIR

module load UHTS/Aligner/bowtie/1.2.0
module load UHTS/Analysis/samtools/1.10

# Mapping to genome
# Files are assumed to be uncompressed
for x in $(ls -d *RNA.fastq); do
bowtie \
 -S \
 -t\
 -p 4 \
 -v 1 \
 -m 1 \
 --best\
 --strata \
 $ADIR \
 $x \
 2> $(basename ${x} .fastq)_GRCh38_log.txt | \
 samtools view \
 -h \
 -F 4 \
 -b > $(basename ${x} .fastq)_GRCh38.bam; done

# Sort the BAM file
for x in $(ls -d *.bam); do
samtools sort \
 -@ 4 \
 $x \
 -o $(basename ${x} .bam)_sorted.bam; done

rm *GRCh38.bam

mkdir $BDIR
mv *.bam $BDIR
