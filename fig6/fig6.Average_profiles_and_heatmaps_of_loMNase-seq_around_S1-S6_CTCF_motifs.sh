#!/usr/bin/env bash

for i in S1 S2 S3 S4 S5 S6; do sed 1d figS2.sameXChIP_differNChIP_${i}.txt | sort -k5gr | cut -f 1-6 >sameXChIP_differNChIP_${i}.bed; done

mkdir -p computeMatrix_outFile/ computeMatrix_outPlot/ computeMatrix_outFileMatrix/

region_file="sameXChIP_differNChIP_S?.bed"

region_label="S1 S2 S3 S4 S5 S6"

samples_label="DMSO TPA"

prefix="K562_loMNase"

bws="K562_loMNase_DMSO.hg38.bw K562_loMNase_TPA.hg38.bw"

computeMatrix reference-point --referencePoint center -S ${bws} -R ${region_file} -b 1000 -a 1000 --missingDataAsZero --averageTypeBins mean --binSize 20 -o computeMatrix_outFile/${prefix}_matrix.gz --outFileNameMatrix computeMatrix_outFileMatrix/${prefix}.matrix -p 80

plotProfile -m computeMatrix_outFile/${prefix}_matrix.gz -out computeMatrix_outPlot/${prefix}_profile.pdf --plotHeight 7 --plotWidth 10 --refPointLabel "center" --colors "#941416" "#C29B39" "#7FC87F" "grey" "#9f79d3" "#386CAF" --samplesLabel ${samples_label} --regionsLabel ${region_label}

plotHeatmap --sortRegions no -m computeMatrix_outFile/${prefix}_matrix.gz -out computeMatrix_outPlot/${prefix}_heatmap.pdf --colorList "white, #08468C" --refPointLabel "center" --heatmapWidth 3 --heatmapHeight 12 --samplesLabel ${samples_label} --whatToShow "heatmap and colorbar" --regionsLabel ${region_label} --boxAroundHeatmaps no --zMax 0.5