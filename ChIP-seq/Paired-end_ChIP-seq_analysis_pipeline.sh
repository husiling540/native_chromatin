#!/usr/bin/env bash

function usage() {
  echo -e "usage : Paired-end_ChIP-seq_analysis_pipeline.sh -i chip_prefixs -p thread -w workdir [-h]"
  #echo -e "for example : Paired-end_ChIP-seq_analysis_pipeline.sh -i "K562_NChIP_CTCF_75mMNaCl_R1 K562_NChIP_CTCF_75mMNaCl_R2 K562_NChIP_CTCF_75mMNaCl_R3" -p 8 -w ~/data"
  echo -e "Use option -h|--help for more information"
}

function help() {
  usage
  echo
  echo "Analyze Paired-end ChIP-seq data. Firstly, you need to create a folder to use as your working path and create a \"1rawdata\" folder inside it to put in your fq.gz file."
  echo "---------------"
  echo "OPTIONS"
  echo
  echo "   -i|--input chip_prefixes : prefixes of ChIP-seq samples"
  echo "   -p|--thread : number of threads during processing"
  echo "   -w|--workdir : workdir"
  echo "   [-h|--help] : help"
  exit
}

if [ $# -lt 1 ]; then
  usage
  exit
fi

while getopts ":i:p:w:h" OPT; do
  case $OPT in

  i) chip_prefixes=$OPTARG ;;
  p) thread=$OPTARG ;;
  w) workdir=$OPTARG ;;
  h) help ;;
  \?)
    echo "Invalid option: -$OPTARG" >&2
    usage
    exit 1
    ;;
  :)
    echo "Option -$OPTARG requires an argument." >&2
    usage
    exit 1
    ;;
  esac
done

# check parameter
if [ -z "${chip_prefixes}" -o -z "${thread}" -o -z "${workdir}" ]; then
  printf "[ERROR] $(date '+%F %T') following parameters is empty:\n-i=${chip_prefixes}\n-p=${thread}\n-w=${workdir}\n"
  exit 1
fi

echo -e "input chip_prefixes: ${chip_prefixes}\nthread: ${thread}\nworkdir: ${workdir}\n============================================"

echo -e "Paired-end ChIP-seq data analysis, start!\n"

cd ${workdir}/

mkdir -p 2fastqc 3trim_adapter 4mapping

for chip_prefix in ${chip_prefixes}; do
  echo -e "${chip_prefix}\t1raw:\t$(zcat 1rawdata/${chip_prefix}_1.fq.gz | awk 'END{print NR/4}')" >${workdir}/Statistics_${chip_prefix}.txt

  echo -e "${chip_prefix} fastqc ======================"
  fastqc -o 2fastqc -q -f fastq 1rawdata/${chip_prefix}*gz -t ${thread}

  echo -e "${chip_prefix} trim_illumina_adapter ======================"
  trim_galore -q 15 --phred33 --stringency 1 --length 10 -e 0.1 --gzip -o 3trim_adapter --paired 1rawdata/${chip_prefix}_1.fq.gz 1rawdata/${chip_prefix}_2.fq.gz
  echo -e "${chip_prefix}\t2trim:\t$(zcat 3trim_adapter/${chip_prefix}_1_val_1.fq.gz | awk 'END{print NR/4}')" >>${workdir}/Statistics_${chip_prefix}.txt
  cutadapt --length 100 -j ${thread} -o 3trim_adapter/${chip_prefix}_trimed_1.fq 3trim_adapter/${chip_prefix}_1_val_1.fq.gz
  cutadapt --length 100 -j ${thread} -o 3trim_adapter/${chip_prefix}_trimed_2.fq 3trim_adapter/${chip_prefix}_2_val_2.fq.gz

  echo -e "${chip_prefix} mapping ======================"
  bowtie2 -x /ssd/index/bowtie2/hg38XX -p ${thread} --no-unal --no-discordant --no-mixed --very-sensitive --score-min L,0,-0.4 -X 1000 -1 3trim_adapter/${chip_prefix}_trimed_1.fq -2 3trim_adapter/${chip_prefix}_trimed_2.fq 2>4mapping/log.${chip_prefix}_bowtie2.txt | samtools view -@ ${thread} -h -q 10 -F 3852 -f 2 | egrep -v chrM | samtools sort -@ ${thread} -o 4mapping/${chip_prefix}.bam

  nQ10=$(samtools view -@ ${thread} -c 4mapping/${chip_prefix}.bam | awk '{print $1/2}')
  echo -e "${chip_prefix}\t3mapping:\t${nQ10}" >>${workdir}/Statistics_${chip_prefix}.txt

  echo -e "${chip_prefix} remove_dup ======================"
  /mnt/disk3/husiling/private/miniconda3/bin/picard MarkDuplicates -REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT --I 4mapping/${chip_prefix}.bam --O 4mapping/${chip_prefix}_rmdup.bam --M 4mapping/${chip_prefix}_rmdup.metrics --OPTICAL_DUPLICATE_PIXEL_DISTANCE 10000 --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 --TAGGING_POLICY All

  index=$(less ./1rawdata/${chip_prefix}_1.fq.gz | head -n 1 | cut -f 2 -d " ")
  READ_PAIR_DUPLICATES=$(cat 4mapping/${chip_prefix}_rmdup.metrics | grep -A 2 "LIBRARY" | cut -f 7-8 | awk 'NR==2{print $1}')
  READ_PAIR_OPTICAL_DUPLICATES=$(cat 4mapping/${chip_prefix}_rmdup.metrics | grep -A 2 "LIBRARY" | cut -f 7-8 | awk 'NR==2{print $2}')
  PCR_dup_rate=$(echo "scale=2; 100*($READ_PAIR_DUPLICATES-$READ_PAIR_OPTICAL_DUPLICATES)/$nQ10" | bc)
  optical_dup_rate=$(echo "scale=2; 100*$READ_PAIR_OPTICAL_DUPLICATES/$nQ10" | bc)

  echo -e "${chip_prefix}\t4rmdup:\t$(samtools view -c 4mapping/${chip_prefix}_rmdup.bam | awk '{print $1/2}')" >>${workdir}/Statistics_${chip_prefix}.txt
  echo -e "${chip_prefix}\t4rmdup_PCR_dup_rate:\t${PCR_dup_rate}%" >>${workdir}/Statistics_${chip_prefix}.txt
  echo -e "${chip_prefix}\t4rmdup_optical_dup_rate:\t${optical_dup_rate}%" >>${workdir}/Statistics_${chip_prefix}.txt

  echo -e "${chip_prefix} convert BAM to BigWig ======================"
  samtools sort -@ ${thread} -n -o 4mapping/tmp.${chip_prefix}_rmdup.sort.bam 4mapping/${chip_prefix}_rmdup.bam
  bamToBed -bedpe -i 4mapping/tmp.${chip_prefix}_rmdup.sort.bam | awk '{print $1,$2,$6}' OFS='\t' | awk '$1!~/chr[CLMT]/ && $2>=0' | sort -k1,1 -k2,2n -S 10G >4mapping/${chip_prefix}_rmdup.bed

  nchr_human=$(wc -l 4mapping/${chip_prefix}_rmdup.bed | cut -f 1 -d " ")
  genomeCoverageBed -bg -i 4mapping/${chip_prefix}_rmdup.bed -g /ssd/genome/hg38_chromsize_real.txt | awk -v nchr_human=${nchr_human} '$4!=0{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/nchr_human}' | awk '{$4/=1;print}' OFS='\t' | sort -k1,1 -k2,2n -S 10G >4mapping/${chip_prefix}.hg38.bgr

  bedGraphToBigWig 4mapping/${chip_prefix}.hg38.bgr /ssd/genome/hg38_chromsize_real.txt 4mapping/${chip_prefix}.hg38.bw

  rm 4mapping/tmp.${chip_prefix}_rmdup.sort.bam
done

cd ${workdir}/