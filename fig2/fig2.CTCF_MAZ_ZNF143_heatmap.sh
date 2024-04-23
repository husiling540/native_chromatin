#!/usr/bin/env bash

mkdir -p computeMatrix_outFile computeMatrix_outFileMatrix computeMatrix_outPlot

#CTCF
prefix="K562_NChIP_XChIP_CTCF"
region_file="CTCF_shared_peak_summits.bed CTCF_XChIP_uniq_peak_summits.bed"
bws="K562_XChIP_CTCF_merge.hg38.bw K562_NChIP_CTCF_150mMNaCl_merge.hg38.bw"
samples_label="XChIP_CTCF NChIP_CTCF"
regionsLabel="CTCF_shared" "XChIP_uniq"

#MAZ
prefix="K562_NChIP_XChIP_MAZ"
region_file="MAZ_shared_peak_summits.bed MAZ_XChIP_uniq_peak_summits.bed MAZ_NChIP_uniq_peak_summits.bed"
bws="K562_ChIP-seq_MAZ_merge.hg38.bw K562_NChIP_MAZ_50mM_merge.hg38.bw"
samples_label="XChIP_MAZ NChIP_MAZ"
regionsLabel="MAZ_shared" "XChIP_uniq" "NChIP_uniq"

#ZNF143
prefix="K562_NChIP_XChIP_ZNF143"
region_file="ZNF143_shared_peak_summits.bed ZNF143_XChIP_uniq_peak_summits.bed ZNF143_NChIP_uniq_peak_summits.bed"
bws="K562_ChIP-seq_ZNF143_Stanford_merge.hg38.bw K562_NChIP_ZNF143_merge.hg38.bw"
samples_label="XChIP_ZNF143 NChIP_ZNF143"
regionsLabel="ZNF143_shared" "XChIP_uniq" "NChIP_uniq"

#plot
computeMatrix reference-point --referencePoint center -S ${bws} -R ${region_file} -b 1000 -a 1000 --missingDataAsZero --averageTypeBins mean --binSize 20 -o computeMatrix_outFile/${prefix}_matrix.gz --outFileNameMatrix computeMatrix_outFileMatrix/${prefix}.matrix -p 80

plotHeatmap -m computeMatrix_outFile/${prefix}_matrix.gz -out computeMatrix_outPlot/${prefix}_heatmap.pdf --colorList "white, #08468C" --refPointLabel "center" --heatmapWidth 3 --heatmapHeight 8 --samplesLabel ${samples_label} --whatToShow "heatmap and colorbar" --regionsLabel ${regionsLabel}
