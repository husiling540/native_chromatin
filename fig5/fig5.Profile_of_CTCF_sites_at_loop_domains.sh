#!/usr/bin/env bash

loop_dir="/mnt/disk3/hsl/private/data/Hi-C/K562/K562_GSE63525/K562_MboI_merge_validPairs_hicpro/hic_results/data/K562_MboI_merge_validPairs/K562_MboI_hiccups_5kb10kb25kb"
native_peak="K562_native_sites.bed"
nonnative_peak="K562_nonnative_sites.bed"
native_count=$(wc -l <${native_peak})
nonnative_count=$(wc -l <${nonnative_peak})

grep -v "#" "$loop_dir"/merged_loops.bedpe | awk '$6-$2>=20000 {print "chr"$1,($2+$3)/2-((($5+$6)/2-($2+$3)/2))/10,($5+$6)/2+((($5+$6)/2-($2+$3)/2))/10,"loop"NR}' OFS='\t' | bedtools makewindows -b - -n 120 -i srcwinnum | sort -k1,1 -k2,2n >K562_loops.120bins.bed

bedtools intersect -a K562_loops.120bins.bed -b ${native_peak} -c | awk -v count=$native_count '{print $1"\t"$2"\t"$3"\t"$4"\t"$5/count*1000}' >K562_loops.120bins.native_sites.bed

bedtools intersect -a K562_loops.120bins.bed -b ${nonnative_peak} -c | awk -v count=$nonnative_count '{print $1"\t"$2"\t"$3"\t"$4"\t"$5/count*1000}' >K562_loops.120bins.nonnative_sites.bed

Rscript fig5.Profile_of_CTCF_sites_at_loop_domains.R