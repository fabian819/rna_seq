#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --partition=pall
#SBATCH --mem-per-cpu=8000M

FDIR="./data/fastq/"
ADIR="/data/users/fgribi/rna_seq/data/annotations/GRCh38_r_t_sno_sn_RNA_ENSEMBL_NCBI_GtRNAdb/GRCh38_r_t_sno_sn_RNA_ENSEMBL_NCBI_GtRNAdb"

cd $FDIR

module load UHTS/Aligner/bowtie/1.2.0

# Mapping to undesired RNAs
# Files are assumed to be uncompressed
for x in $(ls -d *tr.fastq); do
bowtie \
 -S \
 -t\
 -p 4\
 $ADIR \
 $x \
 --un $(basename ${x} .fastq)_no_r_t_sno_sn_RNA.fastq 2> $(basename ${x} .fastq)_no_r_t_sno_sn_RNA_log.txt > /dev/null; done
