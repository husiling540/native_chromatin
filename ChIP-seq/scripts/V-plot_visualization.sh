#!/usr/bin/env bash
##-------------------------
# @Author: Siling Hu
# @Date: 2023/01/01
##-------------------------

mkdir -p ${workdir}/Fig1/vplot
cd ${workdir}/Fig1/vplot

samtools sort -@ 80 -n ${workdir}/Data_preprocessing/4mapping/K562_loMNase_rmdup.bam -o K562_loMNase_sort.bam

bamToBed -bedpe -i K562_loMNase_sort.bam | awk '$1!~/chr[CLMT]/{printf "%s\t%d\t%d\t%d\n",$1,int($2/2+$6/2)-1,int($2/2+$6/2),$6-$2}' |sort -k1,1 -k2,2n -S 16G > K562_loMNase.midpoint.bed &

center=(CTCF ZNF143 MAZ loMNase)
center_file=(/mnt/disk3/husiling/private/data/native_chromatin_summary/Hs_CTCF_motif/Hs_CTCF_conventional.motif /mnt/disk3/husiling/private/data/native_chromatin_summary/K562_ZNF143_motif/K562_ZNF143_motif.bed /mnt/disk3/husiling/private/data/native_chromatin_summary/K562_MAZ_motif/K562_MAZ_HOCOMOCO_motif.bed /mnt/disk3/husiling/private/data/native_chromatin_summary/Data_preprocessing/5macs2/K562_loMNase_macs2/K562_loMNase_peaks.rmBlacklist.narrowPeaksummits.bed)

for i in "${!center[@]}"; do
    echo ${center[i]} ${center_file[i]} 

    awk '{$4=NR;print}' OFS="\t" ${center_file[i]}|cut -f 1-6|bedtools closest -a K562_loMNase.midpoint.bed -b - -D b -t first |awk '!($5=="." && $11=="-1")' > K562_loMNase.midpoint_to_${center[i]}.bed &
done

for sample in CTCF ZNF143 MAZ loMNase; do  Rscript V-plot_visualization.R ${sample} ; done