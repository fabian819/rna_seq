#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --partition=pall
#SBATCH --mem-per-cpu=8000M

ADIR=' ./data/annotations'

module load UHTS/Aligner/bowtie/1.2.0

cd $ADIR

# For the "undesired" RNAs
bowtie-build GRCh38_r-t-sno-sn-RNA.fa GRCh38_r_t_sno_sn_RNA_ENSEMBL_NCBI_GtRNAdb

# For the genome
bowtie-build Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38.dna.primary_assembly

# For the transcriptome
bowtie-build GRCh38_APPRIS_CDS_18.fa GRCh38_APPRIS_CDS_18

# Change APPRIS file to single line
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < GRCh38_APPRIS_CDS_18.fa > GRCh38_APPRIS_CDS_18_SingleLine.fa
