#!/usr/bin/env bash

head -n 17000 K562_conventional_motif.sortbyXChIP.bed >first17000_peak_motif.bed

head -n 34000 K562_conventional_motif.sortbyXChIP.bed | tail -n 17000 >second17000_peak_motif.bed

head -n 51000 K562_conventional_motif.sortbyXChIP.bed | tail -n 17000 >third17000_peak_motif.bed

mkdir -p computeMatrix_outFile computeMatrix_outFileMatrix computeMatrix_outPlot

bws="K562_XChIP_CTCF_merge.hg38.bw K562_NChIP_CTCF_75mMNaCl_merge.hg38.bw K562_NChIP_CTCF_150mMNaCl_merge.hg38.bw K562_NChIP_CTCF_225mMNaCl_merge.hg38.bw "

samples_label="X-ChIP N-ChIP_75mM N-ChIP_150mM N-ChIP_225mM"

prefix=(The_1st_tertile The_2nd_tertile The_3rd_tertile)

region_file=(first17000_peak_motif.bed second17000_peak_motif.bed third17000_peak_motif.bed)

for i in "${!prefix[@]}"; do

  computeMatrix reference-point --referencePoint center -S ${bws} -R ${region_file[i]} -b 1000 -a 1000 --missingDataAsZero --averageTypeBins mean --binSize 20 -o computeMatrix_outFile/${prefix[i]}_matrix.gz --outFileNameMatrix computeMatrix_outFileMatrix/${prefix[i]}.matrix -p 80

  plotProfile -m computeMatrix_outFile/${prefix[i]}_matrix.gz -out computeMatrix_outPlot/${prefix[i]}_profile.pdf --plotHeight 7 --plotWidth 10 --refPointLabel "center" --perGroup --colors "#386CAF" "#7FC87F" "#fc9533" "#9f79d3" "#C09CC2" "#7F7D7E" --samplesLabel ${samples_label} --regionsLabel "merged"

  plotHeatmap -m computeMatrix_outFile/${prefix[i]}_matrix.gz -out computeMatrix_outPlot/${prefix[i]}_heatmap_blue.pdf --colorList "white, #08468C" --refPointLabel "center" --heatmapWidth 2 --heatmapHeight 6.2 --samplesLabel ${samples_label} --whatToShow "heatmap and colorbar" --regionsLabel ${prefix[i]} --sortRegions no --zMax 6 --boxAroundHeatmaps no

done
