#!/usr/bin/env bash

GM12878_dir=/mnt/disk3/hsl/private/data/Hi-C/GM12878/ENCSR410MDC_MboI/GM12878_MboI_merge_validPairs_hicpro/hic_results/data/GM12878_MboI_merge_validPairs/

K562_dir=/mnt/disk3/hsl/private/data/Hi-C/K562/K562_GSE63525/K562_MboI_merge_validPairs_hicpro/hic_results/data/K562_MboI_merge_validPairs/

HeLaS3_dir=/mnt/disk3/hsl/private/data/Hi-C/HeLaS3_HiC_ENCSR693GXU/HeLa-S3_DpnII_merge_validPairs_hicpro/hic_results/data/HeLa-S3_DpnII_merge_validPairs/

IMR90_dir=/mnt/disk3/hsl/private/data/Hi-C/IMR90_HiC_GSE63525/IMR90_MboI_hicpro/hic_results/data/IMR90_MboI/

HUVEC_dir=/mnt/disk3/hsl/private/data/Hi-C/HUVEC_HiC_GSE63525/HUVEC_MboI_hicpro/hic_results/data/HUVEC_MboI/

NHEK_dir=/mnt/disk3/hsl/private/data/Hi-C/NHEK_HiC_GSE63525/NHEK_MboI_hicpro/hic_results/data/NHEK_MboI/

KBM7_dir=/mnt/disk3/hsl/private/data/Hi-C/KBM7_HiC_GSE63525/KBM7_MboI_hicpro/hic_results/data/KBM7_MboI/

HepG2_dir=/mnt/disk3/hsl/private/data/Hi-C/HepG2_GSM3478955/2Analysis/HepG2_MboI_hicpro/hic_results/data/HepG2_MboI/

H1hESC_dir=/mnt/disk3/hsl/private/data/Micro-C/H1hESC_4DNES21D8SP8_validpairs_hicpro/hic_results/data/H1hESC_4DNES21D8SP8_validpairs/

A549_dir=/mnt/disk3/hsl/private/data/Hi-C/A549_HiC_ENCSR662QKG/A549_MboI_merge_validPairs_hicpro/hic_results/data/A549_MboI_merge_validPairs

for cell in GM12878 K562 HeLaS3 IMR90 HUVEC NHEK KBM7 HepG2 A549 H1hESC; do
  tmp=${cell}_dir

  eval dir=$(echo \$$tmp)

  echo ${dir}/TAD_arrowhead_10kb/10000_blocks.bedpe

  grep -v "#" ${dir}/TAD_arrowhead_10kb/10000_blocks.bedpe | awk -v name=$cell '$2>=5000{print "chr"$1"\t"$2-5000"\t"$2+5000"\t"name"\nchr"$1"\t"$3-5000"\t"$3+5000"\t"name}' | sort -k1,1 -k2,2n -k3,3n -u >TAD_10kb_arrowhead_boundaries_${cell}.bed
  wc -l TAD_10kb_arrowhead_boundaries_${cell}.bed
done

for num in $(seq 1 10); do

  cat TAD_10kb_arrowhead_boundaries_*.bed | sort -k1,1 -k2,2n -S 9G | bedtools merge -i - -c 4 -o collapse,count_distinct | awk -v num=$num '$5==num' >TAD_10kb_arrowhead_boundaries.${num}cell.bed

  awk 'length($1)<=5{print $1"\t"int(($2+$3)/2+0.5)-1"\t"int(($2+$3)/2+0.5)"\t"$4"\t"$5}' TAD_10kb_arrowhead_boundaries.${num}cell.bed >TAD_10kb_arrowhead_boundaries.mid.${num}cell.bed

  bedtools slop -i TAD_10kb_arrowhead_boundaries.mid.${num}cell.bed -g /ssd/genome/hg38_chromsize_real.txt -b 500000 | awk -v num=$num '{$4=num"cell";print}' OFS="\t" >TAD_10kb_arrowhead_boundaries.mid_flank500kb.${num}cell.bed

  bedtools makewindows -b TAD_10kb_arrowhead_boundaries.mid_flank500kb.${num}cell.bed -n 100 -i srcwinnum >TAD_10kb_arrowhead_boundaries.mid_flank500kb_100bin.${num}cell.bed

  rm TAD_10kb_arrowhead_boundaries.mid.${num}cell.bed TAD_10kb_arrowhead_boundaries.mid_flank500kb.${num}cell.bed
done

native_peak=K562_native_sites.bed
nonnative_peak=K562_nonnative_sites.bed

native_count=$(wc -l ${native_peak} | cut -f 1 -d " ")
nonnative_count=$(wc -l ${nonnative_peak} | cut -f 1 -d " ")

for num in $(seq 1 10); do
  bedtools intersect -a TAD_10kb_arrowhead_boundaries.mid_flank500kb_100bin.${num}cell.bed -b ${native_peak} -c | awk -v count=$native_count '{print $1"\t"$2"\t"$3"\t"$4"\t"$5/count*1000}' >native_CTCF_sites_at_domain_boundaries_${num}cell.bed
done

for num in $(seq 1 10); do
  bedtools intersect -a TAD_10kb_arrowhead_boundaries.mid_flank500kb_100bin.${num}cell.bed -b ${nonnative_peak} -c | awk -v count=$nonnative_count '{print $1"\t"$2"\t"$3"\t"$4"\t"$5/count*1000}' >nonnative_CTCF_sites_at_domain_boundaries_${num}cell.bed
done

Rscript fig5.Positional_enrichment_of_CTCF_sites_at_domain_boundaries.R