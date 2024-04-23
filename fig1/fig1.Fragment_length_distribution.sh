#!/usr/bin/env bash

#ATAC-seq
samtools view -@ 80 -f 32 native_chromatin_summary/K562_public_data/4mapping/K562_ATAC-seq_PE_rmdup.bam | grep -v chrM | cut -f 9 >readlength_ATAC.txt

#MNase-seq
for i in 5.4U 20.6U 79.2U 304U; do samtools view -@ 80 -f 32 native_chromatin_summary/K562_public_data/4mapping/K562_MNase-seq_${i}_PE_rmdup.bam | cut -f 9 >readlength_MNase_${i}.txt; done

#DNase-seq
samtools view -@ 80 -f 32 native_chromatin_summary/K562_public_data/4mapping/K562_DNase-seq_PE_rmdup.bam | grep -v chrM | cut -f 9 >readlength_DNase.txt

#loMNase-seq
samtools view -@ 80 -f 32 native_chromatin_summary/Data_preprocessing/4mapping/K562_loMNase_rmdup.bam | cut -f 9 >readlength_loMNase.txt

Rscript fig1.Fragment_length_distribution.R