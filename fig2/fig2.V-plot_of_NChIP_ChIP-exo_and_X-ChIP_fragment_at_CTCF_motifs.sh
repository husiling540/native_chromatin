#!/usr/bin/env bash

for sample in K562_NChIP_CTCF_150mMNaCl K562_XChIP_CTCF K562_ChIP-exo_CTCF; do

  bamToBed -bedpe -i ${sample}_merge.sort.bam | awk '{printf "%s\t%d\t%d\t%d\n",$1,int($2/2+$6/2)-1,int($2/2+$6/2),$6-$2}' | sort -k1,1 -k2,2n -S 16G >${sample}.midpoint.bed

  bedtools closest -a ${sample}.midpoint.bed -b Hs_CTCF_conventional.motif -D b -t first | awk '!($5=="." && $11=="-1")' >${sample}.midpoint_to_motif.bed

  Rscript fig2.V-plot_visualization.R ${sample}

done
