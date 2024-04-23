#!/usr/bin/env bash

#K562_ATAC-seq.midpoint
samtools sort -@ 80 -n -o K562_ATAC-seq.sort.bam K562_ATAC-seq_PE_rmdup.bam
bamToBed -bedpe -i K562_ATAC-seq.sort.bam | awk '$1!~/chr[CLMT]/{printf "%s\t%d\t%d\t%d\n",$1,int($2/2+$6/2)-1,int($2/2+$6/2),$6-$2}' | sort -k1,1 -k2,2n --parallel=80 --buffer-size=50% >K562_ATAC.midpoint.bed

#K562_DNase.midpoint
samtools sort -@ 80 -n -o K562_DNase.sort.bam K562_DNase-seq_PE_rmdup.bam
bamToBed -bedpe -i K562_DNase.sort.bam | awk '$1!~/chr[CLMT]/{printf "%s\t%d\t%d\t%d\n",$1,int($2/2+$6/2)-1,int($2/2+$6/2),$6-$2}' | sort -k1,1 -k2,2n --parallel=80 --buffer-size=50% >K562_DNase.midpoint.bed

#K562_loMNase.midpoint
samtools sort -@ 80 -n -o K562_loMNase_sort.bam K562_loMNase_rmdup.bam
bamToBed -bedpe -i K562_loMNase_sort.bam | awk '$1!~/chr[CLMT]/{printf "%s\t%d\t%d\t%d\n",$1,int($2/2+$6/2)-1,int($2/2+$6/2),$6-$2}' | sort -k1,1 -k2,2n -S 16G >K562_loMNase.midpoint.bed

center=(loMNase DNase ATAC)
center_file=(K562_loMNase_peaks.rmBlacklist.narrowPeaksummits.bed K562_DNase-seq_PE_peaks.rmBlacklist.narrowPeaksummits.bed K562_ATAC-seq_PE_peaks.rmBlacklist.narrowPeaksummits.bed)

for i in "${!center[@]}"; do
  echo ${center[i]} ${center_file[i]}
  awk '{$4=NR;print}' OFS="\t" ${center_file[i]} | cut -f 1-6 | bedtools closest -a K562_${center[i]}.midpoint.bed -b - -D b -t first | awk '!($5=="." && $11=="-1")' >K562_${center[i]}.midpoint_to_${center[i]}_summits.bed
done

awk '$11^2<=100' K562_DNase.midpoint_to_DNase_summits.bed >violin_K562_DNase_fragment_length.txt
awk '$11^2<=100' K562_ATAC.midpoint_to_ATAC_summits.bed >violin_K562_ATAC_fragment_length.txt
awk '$11^2<=100' K562_loMNase.midpoint_to_loMNase_summits.bed >violin_K562_loMNase_fragment_length.txt

Rscript fig1.Comparing_the_size_difference_of_fragments_containing_peak_summits.R
