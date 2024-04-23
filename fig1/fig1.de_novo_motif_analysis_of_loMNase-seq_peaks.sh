#!/usr/bin/env bash

flanking=200
count=3000

awk '$1!~/[_,chrM,chrL]/{a=$2;$2=a+$10;$3=a+$10+1;print}' OFS="\t" K562_loMNase_peaks.rmBlacklist.narrowPeak | sort -k1,1 -k2,2n | bedtools intersect -a - -b K562_loMNase.hg38.bgr -wa -wb -sorted | sort -k14gr | head -n ${count} | awk -v flank=${flanking} '$1!~/[_,M,L]/{a=$2;$2=a-flank;$3=a+flank;print}' OFS="\t" >K562_loMNase-seq_peaks.bgrtop${count}_F${flanking}bp.bed

bedtools getfasta -fi hg38XX.fa -bed K562_loMNase-seq_peaks.bgrtop${count}_F${flanking}bp.bed -fo K562_loMNase-seq_peaks.bgrtop${count}_F${flanking}bp.fa

meme-chip -meme-p 80 -order 0 -minw 6 -maxw 25 -meme-nmotifs 5 -oc meme-chip_K562_loMNase-seq_peaks.bgrtop${count}_F${flanking}bp_HOCOMOCO -db /mnt/disk3/hsl/private/reference/Motif_Databases/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme K562_loMNase-seq_peaks.bgrtop${count}_F${flanking}bp.fa

cd meme-chip_K562_loMNase-seq_peaks.bgrtop3000_F200bp_HOCOMOCO/

for i in 1 2 3; do
  count=$(grep -A 3 "MEME-${i} position-specific probability matrix" meme_out/meme.txt | grep 'letter-probability matrix' | grep 'w=' | awk -F 'w=' '{print $2}' | cut -d " " -f 2)
  count2=$(echo -e $count + 2 | bc)
  grep -A ${count2} "MEME-${i} position-specific probability matrix" meme_out/meme.txt >meme_out_MOTIF_${i}.txt
done
