#!/usr/bin/env bash

#The script Codon_occupancy_cal.sh was kindly provided by the Leidel Lab
#and is available at https://github.com/LeidelLab/Codon_occupancy_cal


./Codon_occupancy_cal.sh \
 /data/users/fgribi/rna_seq/data/annotations/GRCh38_APPRIS_CDS_18_SingleLine.fa \
 /data/users/fgribi/rna_seq/data/bam/RPF_KO_Rep1_clpd_tr_no_r_t_sno_sn_RNA_GRCh38_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_KO_Rep1_Codon_occupancy.txt

./Codon_occupancy_cal.sh \
 /data/users/fgribi/rna_seq/data/annotations/GRCh38_APPRIS_CDS_18_SingleLine.fa \
 /data/users/fgribi/rna_seq/data/bam/RPF_KO_Rep2_clpd_tr_no_r_t_sno_sn_RNA_GRCh38_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_KO_Rep2_Codon_occupancy.txt

./Codon_occupancy_cal.sh \
 /data/users/fgribi/rna_seq/data/annotations/GRCh38_APPRIS_CDS_18_SingleLine.fa \
 /data/users/fgribi/rna_seq/data/bam/RPF_WT_Rep1_clpd_tr_no_r_t_sno_sn_RNA_GRCh38_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep1_Codon_occupancy.txt

./Codon_occupancy_cal.sh \
 /data/users/fgribi/rna_seq/data/annotations/GRCh38_APPRIS_CDS_18_SingleLine.fa \
 /data/users/fgribi/rna_seq/data/bam/RPF_WT_Rep2_clpd_tr_no_r_t_sno_sn_RNA_GRCh38_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep2_Codon_occupancy.txt