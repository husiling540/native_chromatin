#!/usr/bin/env bash

# Randomly select midpoints of 10M reads
. ./get_fixed_random.sh
for i in K562_ATAC K562_DNase K562_loMNase; do shuf --random-source=<(get_fixed_random 1836) -n 10000000 ${i}.midpoint.bed | sort -k1,1 -k2,2n --parallel=80 --buffer-size=50% >sample.${i}.midpoint.bed &;done

# Randomly select 90,000 summits from each specified peak file
for peak in /mnt/disk3/hsl/private/data/native_chromatin_summary/Data_preprocessing/5macs2/K562_loMNase_macs2/K562_loMNase_peaks.rmBlacklist.narrowPeaksummits.bed /mnt/disk3/hsl/private/data/native_chromatin_summary/K562_public_data/5macs2/K562_DNase-seq_PE_macs2/K562_DNase-seq_PE_peaks.rmBlacklist.narrowPeaksummits.bed /mnt/disk3/hsl/private/data/native_chromatin_summary/K562_public_data/5macs2/K562_ATAC-seq_PE_macs2/K562_ATAC-seq_PE_peaks.rmBlacklist.narrowPeaksummits.bed; do
  peak_basename=$(basename ${peak})
  shuf --random-source=<(get_fixed_random 1836) -n 90000 ${peak} | sort -k1,1 -k2,2n --parallel=80 --buffer-size=50% >sample.${peak_basename}
done

center=(loMNase DNase ATAC)
center_file=(sample.K562_loMNase_peaks.rmBlacklist.narrowPeaksummits.bed sample.K562_DNase-seq_PE_peaks.rmBlacklist.narrowPeaksummits.bed sample.K562_ATAC-seq_PE_peaks.rmBlacklist.narrowPeaksummits.bed)

for i in "${!center[@]}"; do
  echo ${center[i]} ${center_file[i]}
  awk '{$4=NR;print}' OFS="\t" ${center_file[i]} | cut -f 1-6 | bedtools closest -a sample.K562_${center[i]}.midpoint.bed -b - -D b -t first | awk '!($5=="." && $11=="-1")' >sample.K562_${center[i]}.midpoint_to_${center[i]}_summits.bed &
done

for sample in loMNase DNase ATAC; do
  Rscript fig1.V-plot_of_K562_loMNase-seq_DNase-seq_or_ATAC-seq_fragments_at_their_peak_summits.R ${sample}
done