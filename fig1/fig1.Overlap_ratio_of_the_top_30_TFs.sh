#!/usr/bin/env bash

for peak_file in /mnt/disk3/hsl/private/reference/ENCODE_datasets/ChIP-seq/K562_ENCODE_bw_bed/bed_narrowpeak/TFs/*k.bed; do

  file=K562_loMNase_peaks.rmBlacklist.narrowPeaksummits.bed

  overlap=$(sort -k1,1 -k2,2n ${peak_file} | bedtools intersect -a ${file} -b - -wa -u | wc -l)

  loMNase_peak_count=$(wc -l ${file} | cut -f 1 -d " ")

  TF_peak_count=$(wc -l ${peak_file} | cut -f 1 -d " ")

  peak_basename=$(basename ${peak_file})

  echo hu | awk -v overlap=$overlap -v loMNase_peak_count=$loMNase_peak_count -v peak_basename=$peak_basename -v TF_peak_count=$TF_peak_count '{printf "%s\t%s\t%.4f\n",peak_basename,TF_peak_count,overlap/loMNase_peak_count}' >>loMNase_overlap_ratio.txt
done

awk -F "normal|CRISPR|_" '{print $2"\t"$5}' loMNase_overlap_ratio.txt | sed 's/.narrowPeak.bed//g' | sort -k4gr >loMNase_overlap_ratio.sim.txt