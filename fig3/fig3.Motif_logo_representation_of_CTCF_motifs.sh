#!/usr/bin/env bash

flanking=30

sort -k13gr N1_peak_motif.bed | head -n 100 | awk -v flanking=$flanking '$6=="+"{print $1"\t"$3-flanking+1"\t"$3+flanking+1"\t"$4"\t0\t"$6} $6=="-"{print $1"\t"$2-flanking-1"\t"$2+flanking-1"\t"$4"\t0\t"$6}' | sort -k1,1 -k2,2n >N1_peak_motif_top100.flank${flanking}bp.bed

sort -k12gr N2_peak_motif.bed | head -n 100 | awk -v flanking=$flanking '$6=="+"{print $1"\t"$3-flanking+1"\t"$3+flanking+1"\t"$4"\t0\t"$6} $6=="-"{print $1"\t"$2-flanking-1"\t"$2+flanking-1"\t"$4"\t0\t"$6}' | sort -k1,1 -k2,2n >N2_peak_motif_top100.flank${flanking}bp.bed

sort -k11gr N3_peak_motif.bed | head -n 100 | awk -v flanking=$flanking '$6=="+"{print $1"\t"$3-flanking+1"\t"$3+flanking+1"\t"$4"\t0\t"$6} $6=="-"{print $1"\t"$2-flanking-1"\t"$2+flanking-1"\t"$4"\t0\t"$6}' | sort -k1,1 -k2,2n >N3_peak_motif_top100.flank${flanking}bp.bed

bedtools intersect -a K562_nonnative_sites_motif.bed -b K562_XChIP_CTCF_merge.hg38.bgr -wa -wb -sorted | sort -k14gr -S 9G | head -n 100 | awk -v flanking=$flanking '$6=="+"{print $1"\t"$3-flanking+1"\t"$3+flanking+1"\t"$4"\t0\t"$6} $6=="-"{print $1"\t"$2-flanking-1"\t"$2+flanking-1"\t"$4"\t0\t"$6}' | sort -k1,1 -k2,2n >XChIPspecific_peak_motif_top100.flank${flanking}bp.bed

for group in N1 N2 N3 XChIPspecific; do

  bedtools getfasta -s -fi /ssd/index/bismark/hg38XX/hg38XX.fa -bed ${group}_peak_motif_top100.flank${flanking}bp.bed -fo ${group}_peak_motif_top100.flank${flanking}bp.fa

  awk 'NR%2==0{print}' ${group}_peak_motif_top100.flank${flanking}bp.fa >${group}_peak_motif_top100.flank${flanking}bp.sequence

done

Rscript fig3.Motif_logo_representation_of_CTCF_motifs.R