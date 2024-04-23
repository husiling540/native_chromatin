#!/usr/bin/env bash

if [ -f U2OS_mitosis_peak_fraction.txt ]; then rm U2OS_mitosis_peak_fraction.txt; fi

for fileA in /mnt/disk3/hsl/private/data/ChIP-seq/public_2019_GR_U2OS_RPE1_mitotic_ChIP-seq_CTCF_histone/5macs2/U2OS_ChIP-seq_prometaphase_CTCF_macs2/U2OS_ChIP-seq_prometaphase_CTCF_peaks.rmBlacklist.narrowPeak; do
  for fileB in CTCF_binding_sites/K562_native_sites_summits.bed CTCF_binding_sites/K562_non-native_sites_summits.bed; do

    fileA_basename=$(basename ${fileA})
    fileB_basename=$(basename ${fileB})

    fileA_count=$(wc -l ${fileA} | cut -f 1 -d " ")
    fileB_count=$(wc -l ${fileB} | cut -f 1 -d " ")

    overlap_count=$(bedtools intersect -a ${fileA} -b ${fileB} -wa -u | wc -l)

    echo nchip | awk -v overlap_count=$overlap_count -v fileA_count=$fileA_count -v fileB_count=$fileB_count -v fileA_basename=${fileA_basename} -v fileB_basename=${fileB_basename} '{printf "%s\t%s\t%s\t%s\t%s\t%.3f\n",fileA_basename,fileB_basename,fileA_count,fileB_count,overlap_count,overlap_count/fileA_count}' | sed -e 's/U2OS_ChIP-seq_//g' -e 's/_CTCF_peaks.rmBlacklist.narrowPeak//g' -e 's/K562_//g' -e 's/_sites_summits.bed//g' >>U2OS_mitosis_peak_fraction.txt
  done
done

if [ -f U2OS_interphase_peak_fraction.txt ]; then rm U2OS_interphase_peak_fraction.txt; fi

for fileA in /mnt/disk3/hsl/private/data/ChIP-seq/public_2019_GR_U2OS_RPE1_mitotic_ChIP-seq_CTCF_histone/5macs2/U2OS_ChIP-seq_interphase_CTCF_macs2/U2OS_ChIP-seq_interphase_CTCF_peaks.rmBlacklist.narrowPeak; do
  for fileB in CTCF_binding_sites/K562_native_sites_summits.bed CTCF_binding_sites/K562_non-native_sites_summits.bed; do

    fileA_basename=$(basename ${fileA})
    fileB_basename=$(basename ${fileB})

    fileA_count=$(wc -l ${fileA} | cut -f 1 -d " ")
    fileB_count=$(wc -l ${fileB} | cut -f 1 -d " ")

    overlap_count=$(bedtools intersect -a ${fileA} -b ${fileB} -wa -u | wc -l)

    echo nchip | awk -v overlap_count=$overlap_count -v fileA_count=$fileA_count -v fileB_count=$fileB_count -v fileA_basename=${fileA_basename} -v fileB_basename=${fileB_basename} '{printf "%s\t%s\t%s\t%s\t%s\t%.3f\n",fileA_basename,fileB_basename,fileA_count,fileB_count,overlap_count,overlap_count/fileA_count}' | sed -e 's/U2OS_ChIP-seq_//g' -e 's/_CTCF_peaks.rmBlacklist.narrowPeak//g' -e 's/K562_//g' -e 's/_sites_summits.bed//g' >>U2OS_interphase_peak_fraction.txt
  done
done

Rscript fig5.Percentage_of_CTCF_sites-Interphase_and_prometaphase.R