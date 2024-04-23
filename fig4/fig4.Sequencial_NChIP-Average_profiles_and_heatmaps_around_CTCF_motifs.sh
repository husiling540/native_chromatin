#!/usr/bin/env bash

mkdir -p computeMatrix_outFile computeMatrix_outFileMatrix computeMatrix_outPlot

computeMatrix reference-point --referencePoint center -S K562_NChIP_CTCF_75mMNaCl_merge.hg38.bw K562_NChIP_CTCF_150mMNaCl_merge.hg38.bw K562_NChIP_CTCF_225mMNaCl_merge.hg38.bw -R K562_CTCF_SUMO_sites.CTCF_motif.bed K562_CTCF_non-SUMO_sites.CTCF_motif.bed -b 500 -a 500 --missingDataAsZero --averageTypeBins mean --binSize 5 -o computeMatrix_outFile/K562_NChIP-3salt_matrix.gz --outFileNameMatrix computeMatrix_outFileMatrix/K562_NChIP_CTCF_different_NaCl.matrix -p 40

plotProfile -m computeMatrix_outFile/K562_NChIP_CTCF_different_NaCl_matrix.gz -out computeMatrix_outPlot/K562_NChIP_CTCF_different_NaCl_profile.pdf --plotHeight 7 --plotWidth 10 --colors "#C36518" "#2A5CA3" --samplesLabel 75mMNaCl 150mMNaCl 225mMNaCl

plotHeatmap -m computeMatrix_outFile/K562_NChIP_CTCF_different_NaCl_matrix.gz -out computeMatrix_outPlot/K562_NChIP_CTCF_different_NaCl.pdf --colorList "white, #08468C" --refPointLabel "" --regionsLabel SUMO non-SUMO --samplesLabel 75mMNaCl 150mMNaCl 225mMNaCl --heatmapWidth 3 --heatmapHeight 8 --xAxisLabel "" --whatToShow "plot and heatmap" --whatToShow "heatmap and colorbar"