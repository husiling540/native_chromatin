#!/usr/bin/env bash

samtools sort -@ 80 -n K562_loMNase_rmdup.bam -o K562_loMNase_sort.bam
bamToBed -bedpe -i K562_loMNase_sort.bam | awk '$1!~/chr[CLMT]/{printf "%s\t%d\t%d\t%d\n",$1,int($2/2+$6/2)-1,int($2/2+$6/2),$6-$2}' | sort -k1,1 -k2,2n -S 16G >K562_loMNase.midpoint.bed &

center=(CTCF ZNF143_MOTIF1 ZNF143_MOTIF2 MAZ)
center_file=(Hs_CTCF_conventional.motif K562_ZNF143_MOTIF1.bed K562_ZNF143_MOTIF2.bed K562_MAZ_HOCOMOCO_motif.bed)

for i in "${!center[@]}"; do
  echo ${center[i]} ${center_file[i]}
  awk '{$4=NR;print}' OFS="\t" ${center_file[i]} | cut -f 1-6 | bedtools closest -a K562_loMNase.midpoint.bed -b - -D b -t first | awk '!($5=="." && $11=="-1")' >K562_loMNase.midpoint_to_${center[i]}.bed &
done

for sample in CTCF ZNF143_MOTIF1 ZNF143_MOTIF2 MAZ; do Rscript fig1.V-plot_visualization.R ${sample}; done