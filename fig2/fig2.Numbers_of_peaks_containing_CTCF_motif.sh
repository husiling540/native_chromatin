#!/usr/bin/env bash

for sample in K562_CUT-Tag_CTCF K562_CUT-RUN_CTCF K562_NChIP-public_CTCF K562_NChIP_CTCF_150mMNaCl; do

  samtools view -@ 80 -h -f 64 -b ${sample}_merge.sort.bam -o ${sample}_merge_R1.bam
  samtools view -@ 80 -h -f 128 -b ${sample}_merge.sort.bam -o ${sample}_merge_R2.bam

done

for i in K562_CUT-Tag_CTCF_merge_R1.bam K562_CUT-Tag_CTCF_merge_R2.bam K562_CUT-RUN_CTCF_merge_R1.bam K562_CUT-RUN_CTCF_merge_R2.bam K562_NChIP-public_CTCF_merge_R1.bam K562_NChIP-public_CTCF_merge_R2.bam K562_NChIP_CTCF_150mMNaCl_merge_R1.bam K562_NChIP_CTCF_150mMNaCl_merge_R2.bam; do

  /home/hsl/softwares/bbmap/reformat.sh in=${i} out=subsampled.${i} samplereads=2100000 sampleseed=1836 overwrite=true

done

for sample in K562_CUT-Tag_CTCF K562_CUT-RUN_CTCF K562_NChIP-public_CTCF K562_NChIP_CTCF_150mMNaCl; do

  samtools view subsampled.${sample}_merge_R1.bam | split -d -l 700000 - subsampled.${sample}_merge_R1_split
  samtools view subsampled.${sample}_merge_R2.bam | split -d -l 700000 - subsampled.${sample}_merge_R2_split

  rm subsampled.${sample}_merge_R1.bam subsampled.${sample}_merge_R2.bam

  for read in R1 R2; do

    samtools view -H subsampled.${sample}_merge_${read}.bam >subsampled.${sample}_merge_${read}.bam.header

    cat subsampled.${sample}_merge_${read}.bam.header subsampled.${sample}_merge_${read}_split00 | samtools view -bh -@ 48 - >subsampled.${sample}_merge_${read}_split_1.bam

    cat subsampled.${sample}_merge_${read}.bam.header subsampled.${sample}_merge_${read}_split00 subsampled.${sample}_merge_${read}_split01 | samtools view -bh -@ 48 - >subsampled.${sample}_merge_${read}_split_2.bam

    cat subsampled.${sample}_merge_${read}.bam.header subsampled.${sample}_merge_${read}_split00 subsampled.${sample}_merge_${read}_split01 subsampled.${sample}_merge_${read}_split02 | samtools view -bh -@ 48 - >subsampled.${sample}_merge_${read}_split_3.bam

    rm subsampled.${sample}_merge_${read}.bam.header subsampled.${sample}_merge_${read}_split00 subsampled.${sample}_merge_${read}_split01 subsampled.${sample}_merge_${read}_split02

  done

  for i in 1 2 3; do

    samtools merge -@ 80 -f subsampled.${sample}_merge_split_${i}.bam subsampled.${sample}_merge_R1_split_${i}.bam subsampled.${sample}_merge_R2_split_${i}.bam

    samtools sort -@ 80 -n -o subsampled.${sample}_merge_split_${i}.sort.bam subsampled.${sample}_merge_split_${i}.bam

    rm subsampled.${sample}_merge_split_${i}.bam subsampled.${sample}_merge_R1_split_${i}.bam subsampled.${sample}_merge_R2_split_${i}.bam

  done
done

#K562_NChIP_CTCF_150mMNaCl
for sample in K562_NChIP_CTCF_150mMNaCl_merge_split_1 K562_NChIP_CTCF_150mMNaCl_merge_split_2 K562_NChIP_CTCF_150mMNaCl_merge_split_3; do
  macs2 callpeak -t subsampled.${sample}.sort.bam -f BAM --keep-dup all -n ${sample} --nomodel -q 0.05 -g hs --outdir macs2/${sample}_macs2/ >macs2/${sample}_macs2.log 2>&1 &
done

#K562_CUT-RUN_CTCF
for sample in K562_CUT-RUN_CTCF_merge_split_1 K562_CUT-RUN_CTCF_merge_split_2 K562_CUT-RUN_CTCF_merge_split_3; do
  macs2 callpeak -t subsampled.${sample}.sort.bam -f BAM --keep-dup all -n ${sample} --nomodel -q 0.05 -g hs --outdir macs2/${sample}_macs2/ >macs2/${sample}_macs2.log 2>&1 &
done

#K562_CUT-Tag_CTCF
for sample in K562_CUT-Tag_CTCF_merge_split_1 K562_CUT-Tag_CTCF_merge_split_2 K562_CUT-Tag_CTCF_merge_split_3; do
  macs2 callpeak -t subsampled.${sample}.sort.bam -f BAM --keep-dup all -n ${sample} --nomodel -q 0.05 -g hs --outdir macs2/${sample}_macs2/ >macs2/${sample}_macs2.log 2>&1 &
done

#K562_NChIP-public_CTCF
for sample in K562_NChIP-public_CTCF_merge_split_1 K562_NChIP-public_CTCF_merge_split_2 K562_NChIP-public_CTCF_merge_split_3; do
  macs2 callpeak -t subsampled.${sample}.sort.bam -f BAM --keep-dup all -n ${sample} --nomodel -q 0.05 -g hs --outdir macs2/${sample}_macs2/ >macs2/${sample}_macs2.log 2>&1 &
done

cd macs2

for sample in K562_CUT-Tag_CTCF_merge_split_1 K562_CUT-Tag_CTCF_merge_split_2 K562_CUT-Tag_CTCF_merge_split_3 K562_CUT-RUN_CTCF_merge_split_1 K562_CUT-RUN_CTCF_merge_split_2 K562_CUT-RUN_CTCF_merge_split_3 K562_NChIP-public_CTCF_merge_split_1 K562_NChIP-public_CTCF_merge_split_2 K562_NChIP-public_CTCF_merge_split_3 K562_NChIP_CTCF_150mMNaCl_merge_split_1 K562_NChIP_CTCF_150mMNaCl_merge_split_2 K562_NChIP_CTCF_150mMNaCl_merge_split_3; do
  sort -k1,1 -k2,2n ${sample}_macs2/${sample}_peaks.narrowPeak | bedtools intersect -a - -b GRCh38_blacklisted_regions_ENCFF356LFX.bed -v >${sample}_macs2/${sample}_peaks.rmBlacklist.narrowPeak
  awk '{s=$2;d=$10;$2=s+d;$3=s+d+1;print}' OFS="\t" ${sample}_macs2/${sample}_peaks.rmBlacklist.narrowPeak >${sample}_macs2/${sample}_peaks.rmBlacklist.narrowPeaksummits.bed
done

if [ -f Ratio_of_peaks_containing_CTCF_motif.txt ]; then rm Ratio_of_peaks_containing_CTCF_motif.txt; fi

for file in *macs2/*.rmBlacklist.narrowPeak; do
  motif=$(bedtools intersect -a ${file} -b ${motif_file} -wa -u | wc -l)
  all=$(wc -l ${file} | cut -f 1 -d " ")
  file_basename=$(basename ${file})
  echo hu | awk -v motif=$motif -v all=$all -v file_basename=$file_basename '{printf "%s\t%s\t%s\t%.4f\n",file_basename,motif,all-motif,motif/all}' | sed -e 's/K562_//g' -e 's/CTCF_merge_split_//g' -e 's/_peaks.rmBlacklist.narrowPeak//g' -e 's/_CTCF_150mMNaCl_merge_split//g' -e 's/_/\t/g' | grep -v XChIP >>Number_of_peaks_containing_CTCF_motif.txt
done

Rscript fig2.Numbers_of_peaks_containing_CTCF_motif.R