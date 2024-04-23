#!/usr/bin/env bash

prefix=K562_native_sites_6library

awk '$5>=6' OFS="\t" K562_native_sites.bed > ${prefix}.bed;done

bedtools intersect -a ${prefix}.bed -b Hs_CTCF_conventional_motif.txt -wa -v|awk '{mid=int($2/2+$3/2);$2=int(mid-75);$3=(mid+75);print}' OFS="\t"> ${prefix}.noCTCFmotif.bed
	
bedtools getfasta -fi /ssd/index/bismark/hg38XX/hg38XX.fa -bed ${prefix}.noCTCFmotif.bed -fo ${prefix}.noCTCFmotif.fa

#Identification of ​V-motif
meme ${prefix}.noCTCFmotif.fa -revcomp -dna -nmotifs 3 -minw 8 -maxw 19 -mod zoops -nostatus -oc meme_${prefix}_noCTCFmotif


count=`grep -A 3 "MEME-1 position-specific probability matrix" meme_K562_native_sites_6library_noCTCFmotif/meme.txt |  grep 'letter-probability matrix'|grep 'w=' | awk -F 'w=' '{print $2}'|cut -d " " -f 2`

count2=`echo -e $count + 2|bc`
	
grep -A ${count2} "MEME-1 position-specific probability matrix" meme_K562_native_sites_6library_noCTCFmotif/meme.txt > meme_K562_native_sites_6library_noCTCFmotif_MOTIF_1.txt

Rscript fig3.CTCF_V-motif.R