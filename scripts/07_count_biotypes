#!/usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --partition=pall
#SBATCH --mem-per-cpu=1000M

module load UHTS/Analysis/subread/2.0.1

BDIR="/data/users/fgribi/rna_seq/data/bam"
ADIR="/data/users/fgribi/rna_seq/data/annotations/Homo_sapiens.GRCh38.108.gtf.gz"

cd $BDIR


# Extract reads mapped to different biotypes

featureCounts \
 -T 8 \
 -t exon \
 -g gene_biotype \
 -a $ADIR \
 -o biotype_counts_rawfile.txt *_GRCh38_sorted.bam

# Extract Biotype and Sample columns
cut -f 1,7-10 biotype_counts_rawfile.txt > biotype_counts_processed.txt