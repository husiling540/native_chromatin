#!/usr/bin/env bash

for i in N1 N2 N3 XChIP; do sed 1d figEV4.sameXChIP_differNChIP_${i}.txt | sort -k5gr >sameXChIP_differNChIP_${i}.bed; done

mkdir -p computeMatrix_outFile/ computeMatrix_outPlot/heatmap computeMatrix_outFileMatrix/

bed_file="sameXChIP_differNChIP_N1.bed sameXChIP_differNChIP_N2.bed sameXChIP_differNChIP_N3.bed sameXChIP_differNChIP_XChIP.bed"
region_label="N1 N2 N3 X-ChIP"
flank=1000
binSize=10
line_colors=("#386CAF" "#7FC87F" "#fc9533" "#9f79d3" "#81CFD9" "#F6C83D")
thread=40
plotHeight=7
plotWidth=10
refPointLabel="motif center"

ls /mnt/disk3/hsl/private/data/native_chromatin_summary/K562_public_data/4mapping/K562_ChIP-seq_*bw | sed -e 's/\///g' -e 's/.bw//g' | sed -e 's/mntdisk3hslprivatedatanative_chromatin_summaryK562_public_data4mapping//g' >bw_K562_public_data_to_plot_profile.txt

for sample in $(cat bw_K562_public_data_to_plot_profile.txt); do
  echo "${sample} /mnt/disk3/hsl/private/data/native_chromatin_summary/K562_public_data/4mapping/${sample}.bw"

  computeMatrix reference-point --referencePoint center -S /mnt/disk3/hsl/private/data/native_chromatin_summary/K562_public_data/4mapping/${sample}.bw -R ${bed_file} -a ${flank} -b ${flank} --missingDataAsZero --averageTypeBins mean --binSize ${binSize} -o computeMatrix_outFile/K562_Xu_TF_${sample}_matrix.gz --outFileNameMatrix computeMatrix_outFileMatrix/K562_Xu_TF_${sample}.matrix -p ${thread}

  plotProfile -m computeMatrix_outFile/K562_Xu_TF_${sample}_matrix.gz -out computeMatrix_outPlot/K562_Xu_TF_${sample}.pdf --regionsLabel ${region_label} --colors ${line_colors[@]} --plotHeight ${plotHeight} --plotWidth ${plotWidth} --refPointLabel "${refPointLabel}"
done

sed -e 's/K562_ChIP-seq_//g' -e 's/.hg38//g' bw_K562_public_data_to_plot_profile.txt >bw_K562_public_data_to_plot_profile.sim.txt

Rscript fig4.Average_profiles.R