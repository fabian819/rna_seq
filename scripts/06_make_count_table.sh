#!/usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --partition=pall
#SBATCH --mem-per-cpu=1000M

module load UHTS/Analysis/subread/2.0.1

BDIR="/data/users/fgribi/rna_seq/data/bam"
ADIR="/data/users/fgribi/rna_seq/data/annotations/Homo_sapiens.GRCh38.108.gtf.gz"

cd $BDIR

# Count reads on CDS
featureCounts \
 -T 8 \
 -t CDS \
 -g gene_id \
 -a $ADIR \
 -o CDS_counts_rawfile.txt *_GRCh38_sorted.bam

# Extract GeneID and Sample columns
cut -f 1,7-10 CDS_counts_rawfile.txt > CDS_counts_processed.txt