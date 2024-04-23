#!/usr/bin/env bash

for prefix in K562_XChIP_CTCF_merge K562_NChIP_CTCF_150mMNaCl_merge K562_CUT-RUN_CTCF_merge K562_CUT-Tag_CTCF_merge K562_ChIP-exo_CTCF_merge; do

  if [ ! -f ${prefix}.sort.bam ]; then samtools sort -@ 80 -n -o ${prefix}.sort.bam ${prefix}.bam; fi

  bamToBed -bedpe -i ${prefix}.sort.bam | awk '{print $1"\t"$2"\t"$6"\t"$7}' | sort -k1,1 -k2,2n -S 16G >${prefix}.bed

  awk '$3-$2 <= 120' ${prefix}.bed | bedtools closest -a - -b /mnt/disk3/hsl/private/data/native_chromatin_summary/Hs_CTCF_motif/Hs_CTCF_conventional.motif -D b -t first | awk '$10=="+"{print $1"\t"$2"\t"$2+1"\t"$4"\tup\n"$1"\t"$3-1"\t"$3"\t"$4"\tdown"}$10=="-"{print $1"\t"$2"\t"$2+1"\t"$4"\tdown\n"$1"\t"$3-1"\t"$3"\t"$4"\tup"}' >${prefix}.border.bed

  grep up ${prefix}.border.bed | sort -k1,1 -k2,2n --parallel=40 --buffer-size=50% >${prefix}.5border.bed

  grep down ${prefix}.border.bed | sort -k1,1 -k2,2n --parallel=40 --buffer-size=50% >${prefix}.3border.bed

  n_fragment=$(wc -l ${prefix}.5border.bed | cut -f 1 -d " ")

  for end in 5 3; do
    genomeCoverageBed -bg -i ${prefix}.${end}border.bed -g /ssd/genome/hg38_chromsize_real.txt | awk -v n_fragment=${n_fragment} '$4!=0{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/n_fragment}' | awk '{$4/=1;print}' OFS='\t' | sort -k1,1 -k2,2n --parallel=40 --buffer-size=50% >${prefix}.${end}border.bgr &&
      bedGraphToBigWig ${prefix}.${end}border.bgr /ssd/genome/hg38_chromsize_real.txt ${prefix}.${end}border.bw
  done
done

bws="K562_XChIP_CTCF_merge.3border.bw K562_XChIP_CTCF_merge.5border.bw K562_NChIP_CTCF_150mMNaCl_merge.3border.bw K562_NChIP_CTCF_150mMNaCl_merge.5border.bw K562_CUT-RUN_CTCF_merge.3border.bw K562_CUT-RUN_CTCF_merge.5border.bw K562_CUT-Tag_CTCF_merge.3border.bw K562_CUT-Tag_CTCF_merge.5border.bw K562_ChIP-exo_CTCF_merge.3border.bw K562_ChIP-exo_CTCF_merge.5border.bw "

samples_label=$(echo $bws | sed -e 's/border.bw//g' -e 's/K562_//g')

computeMatrix reference-point --referencePoint center -S ${bws} -R /mnt/disk3/hsl/private/data/native_chromatin_summary/Hs_CTCF_motif/Hs_CTCF_conventional.motif -b 100 -a 100 --missingDataAsZero --averageTypeBins mean --binSize 1 -o motifcenter.K562_footprint.gz --outFileNameMatrix motifcenter.K562_footprint.Matrix -p 80

plotProfile -m motifcenter.K562_footprint.gz -out motifcenter.K562_footprint.pdf --perGroup --refPointLabel "motif center" --samplesLabel ${samples_label} --regionsLabel "native footprint" --plotHeight 7 --plotWidth 10

Rscript fig2.Scaled_pileup_of_fragment_ends.R