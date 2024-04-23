#!/usr/bin/env bash

for file in K562_CTCF_SUMO_sites.bed K562_CTCF_non-SUMO_sites.bed; do
  echo ${file}
  cut -f 1-3 ${file} | bedtools intersect -a - -b K562_cmotif.salt_tolerance_score.bed -wa -wb | awk '{$13=(int($6-($2+$3)/2)^2)^0.5;print}' OFS="\t" | sort -k1,1 -k2,2n -k13g -S 16G | cut -f 1-12 | awk '!a[$1$2$3]++' >${file%.bed}.salt_tolerance_score.txt
done

Rscript fig4.Sequencial_NChIP-Boxplot_of_S-score.R